import requests
# import dncs
import time
import zipfile
import os
from io import BytesIO

class SimulationClient:
    def __init__(self, base_url="http://localhost:3000", sequence="", samples=100, force_field="amberfb15.xml", grid=4):
        """Initialize the simulation client with form data and API base URL."""
        self.base_url = base_url.rstrip("/")
        self.form_data = {
            "sequence": sequence,
            "samples": samples,
            "force_field": force_field,
            "grid": grid,
        }
        self.process_id = None
        self.status_message = ""
        self.simulation_progress = 0
        self.headers = {"Content-Type": "application/json"}

    def validate_inputs(self):
        """Validate form data before sending requests."""
        if not self.form_data["sequence"]:
            print("Warning: Sequence is empty. Ensure this is valid for the API.")
        if not (1 <= self.form_data["samples"] <= 1000):
            raise ValueError("Samples must be between 1 and 1000.")
        return True

    def send_simulation_request(self):
        """Send POST request to /api/sample to start the simulation."""
        url = f"{self.base_url}/api/sample"
        payload = {
            "sequence": self.form_data["sequence"],
            "sample": self.form_data["samples"],
            "forcefield": self.form_data["force_field"],
            "grid": self.form_data["grid"],
        }
        
        try:
            response = requests.post(url, json=payload, headers=self.headers)
            print(f"Simulation Request HTTP Status Code: {response.status_code}")
            print(f"Simulation Request Raw Response: '{response.text}'")
            print(f"Response Headers: {response.headers}")
            response.raise_for_status()
            
            # Clean and validate the process ID
            self.process_id = response.text.strip()
            if not self.process_id:
                print("❌ Warning: Empty process ID received from server")
                return None
            
            print(f"✅ Process ID received: '{self.process_id}'")
            print(f"   Length: {len(self.process_id)} characters")
            print(f"   Contains only alphanumeric: {self.process_id.isalnum()}")
            
            return self.process_id
        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
            print(f"Response: {response.text}")
        except requests.exceptions.ConnectionError:
            print("Error: Could not connect to server. Is it running?")
        except requests.exceptions.RequestException as err:
            print(f"Error sending request: {err}")
        return None

    def poll_status(self, max_duration=300, poll_interval=3, max_invalid_responses=10):
        """Poll /api/status until completion or timeout, with improved debugging."""
        if not self.process_id:
            raise ValueError("No process ID available. Run send_simulation_request first.")
        
        start_time = time.time()
        invalid_response_count = 0
        
        # Based on test results, we know the correct format
        base_url = f"{self.base_url}/api/status?id={self.process_id}"
        
        print(f"🔄 Starting status polling for process ID: '{self.process_id}'")
        print(f"📍 Using URL: {base_url}")
        
        # Wait a moment before first poll to let the server register the process
        print("⏳ Waiting 3 seconds for server to register the process...")
        time.sleep(3)
        
        while time.time() - start_time < max_duration:
            try:
                # Add headers that might be expected by the server
                headers = {
                    "Accept": "text/plain, application/json, */*",
                    "User-Agent": "Python-SimulationClient/1.0",
                    "Cache-Control": "no-cache",
                    "Pragma": "no-cache"
                }
                
                response = requests.get(base_url, headers=headers, timeout=10)
                
                print(f"\n⏱️ Poll at {time.time() - start_time:.1f}s:")
                print(f"   HTTP Status: {response.status_code}")
                print(f"   Response: '{response.text}'")
                print(f"   Content-Type: {response.headers.get('Content-Type', 'Not specified')}")
                
                # Handle different response codes
                if response.status_code == 200:
                    new_status = response.text.strip() if response.text else ""
                    
                    if new_status and new_status.lower() not in {"", "null", "undefined", "not found html"}:
                        # Reset invalid response counter on success
                        invalid_response_count = 0
                        
                        # Only update if status actually changed
                        if new_status != self.status_message:
                            self.status_message = new_status
                            print(f"✅ Status Updated: '{self.status_message}'")
                            
                            # Update progress based on status (case-insensitive)
                            status_lower = self.status_message.lower()
                            if "waiting" in status_lower or "queue" in status_lower:
                                self.simulation_progress = 10
                            elif "generating" in status_lower or "sample" in status_lower:
                                self.simulation_progress = 33
                            elif "minimizing" in status_lower or "minimize" in status_lower:
                                self.simulation_progress = 66
                            elif "completed" in status_lower or "complete" in status_lower or "done" in status_lower:
                                self.simulation_progress = 100
                                print("🎉 Simulation completed successfully!")
                                return True
                            elif "error" in status_lower or "failed" in status_lower or "fail" in status_lower:
                                print(f"❌ Simulation failed: {self.status_message}")
                                return False
                            
                            print(f"📊 Progress: {self.simulation_progress}%")
                        else:
                            print(f"🔄 Status unchanged: '{self.status_message}'")
                    else:
                        invalid_response_count += 1
                        print(f"⚠️ Empty/invalid response #{invalid_response_count}")
                        
                elif response.status_code == 404:
                    print(f"❌ Process ID '{self.process_id}' not found on server")
                    print("   This could mean:")
                    print("   - The process hasn't been registered yet (will retry)")
                    print("   - The process ID is incorrect")
                    print("   - The process has expired")
                    invalid_response_count += 1
                    
                elif response.status_code == 400:
                    print(f"❌ Bad request: {response.text}")
                    print("   Check if the process ID format is correct")
                    return False
                    
                else:
                    print(f"❌ Unexpected HTTP status: {response.status_code}")
                    print(f"   Response: {response.text}")
                    invalid_response_count += 1
                
            except requests.exceptions.Timeout:
                print("⏱️ Request timeout")
                invalid_response_count += 1
            except requests.exceptions.ConnectionError:
                print("🔌 Connection error")
                invalid_response_count += 1
            except requests.exceptions.RequestException as err:
                print(f"❌ Request error: {err}")
                invalid_response_count += 1
            
            # Check if we've had too many invalid responses
            if invalid_response_count >= max_invalid_responses:
                print(f"❌ Too many invalid responses ({invalid_response_count}). Aborting polling.")
                return False
            
            print(f"⏳ Waiting {poll_interval} seconds before next poll...")
            time.sleep(poll_interval)
        
        print(f"⏱️ Polling timed out after {max_duration} seconds.")
        print(f"Last known status: '{self.status_message}'")
        return False

    def test_status_endpoint(self):
        """Test different status endpoint formats to see which one works."""
        if not self.process_id:
            print("❌ No process ID available. Run send_simulation_request first.")
            return
        
        print(f"🔍 Testing status endpoints for process ID: {self.process_id}")
        
        test_urls = [
            f"{self.base_url}/api/status?id={self.process_id}",
            f"{self.base_url}/api/status/{self.process_id}",
            f"{self.base_url}/api/status?processId={self.process_id}",
            f"{self.base_url}/api/status?pid={self.process_id}",
            f"{self.base_url}/status?id={self.process_id}",
            f"{self.base_url}/status/{self.process_id}"
        ]
        
        for url in test_urls:
            try:
                response = requests.get(url, timeout=5)
                print(f"\n📍 URL: {url}")
                print(f"   Status: {response.status_code}")
                print(f"   Response: '{response.text[:100]}{'...' if len(response.text) > 100 else ''}'")
                print(f"   Content-Type: {response.headers.get('Content-Type', 'Not specified')}")
            except Exception as e:
                print(f"\n📍 URL: {url}")
                print(f"   Error: {e}")

    def load_first_sample(self):
        """Fetch the first sample PDB file and save it."""
        if not self.process_id:
            raise ValueError("No process ID available. Run send_simulation_request first.")
        
        url = f"{self.base_url}/api/pdb?path=Result/{self.process_id}/sample/sample_0000.pdb"
        headers = {
            "Cache-Control": "no-cache",
            "Pragma": "no-cache"
        }
        
        try:
            response = requests.get(url, headers=headers)
            print(f"First Sample HTTP Status Code: {response.status_code}")
            print(f"First Sample Raw Response: {response.text[:100]}...")
            # print(f"First Sample Response Headers: {response.headers}")
            response.raise_for_status()
            with open("sample_0000.pdb", "w") as f:
                f.write(response.text)
            print("Saved first sample to sample_0000.pdb")
            return response.text
        except requests.exceptions.HTTPError as http_err:
            print(f"HTTP error fetching first sample: {http_err}")
            print(f"Response: {response.text}")
        except requests.exceptions.RequestException as err:
            print(f"Error fetching first sample: {err}")
        return None

    def download_data(self):
        """Download all sample files into a ZIP."""
        if not self.process_id:
            raise ValueError("No process ID available. Run send_simulation_request first.")
        
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for i in range(self.form_data["samples"]):
                padded_number = f"{i:04d}"
                sample_url = f"{self.base_url}/api/pdb?path=Result/{self.process_id}/sample/sample_{padded_number}.pdb"
                sample_out_url = f"{self.base_url}/api/pdb?path=Result/{self.process_id}/sample/sampled.out"
                
                try:
                    sample_response = requests.get(sample_url)
                    sample_out_response = requests.get(sample_out_url)
                    
                    print(f"Sample {padded_number} HTTP Status: {sample_response.status_code}")
                    print(f"Sampled.out HTTP Status: {sample_out_response.status_code}")
                    # print(f"Sample Response Headers: {sample_response.headers}")
                    
                    sample_response.raise_for_status()
                    sample_out_response.raise_for_status()
                    
                    if sample_response.text:
                        zip_file.writestr(
                            f"{self.process_id}/sample/sample_{padded_number}.pdb",
                            sample_response.text
                        )
                    if sample_out_response.text:
                        zip_file.writestr(
                            f"{self.process_id}/sample/sampled.out",
                            sample_out_response.text
                        )
                except requests.exceptions.RequestException as err:
                    print(f"Error fetching sample {padded_number}: {err}")
        
        zip_buffer.seek(0)
        zip_path = f"{self.process_id}.zip"
        with open(zip_path, "wb") as f:
            f.write(zip_buffer.read())
        print(f"Saved ZIP file to {zip_path}")

    def run(self):
        """Execute the full simulation workflow."""
        try:
            self.validate_inputs()
            
            print("🚀 Starting simulation...")
            if not self.send_simulation_request():
                print("❌ Failed to start simulation.")
                return False
            
            print(f"✅ Simulation started with process ID: '{self.process_id}'")
            
            # Optional: Test status endpoints if you want to debug
            # print("\n🔍 Testing status endpoints...")
            # self.test_status_endpoint()
            
            print("\n📊 Polling simulation status...")
            if not self.poll_status():
                print("❌ Failed to complete simulation or polling timed out.")
                print("💡 The simulation might still be running on the server.")
                print("   Check the web interface to see the current status.")
                return False
            
            print("\n📁 Fetching first sample...")
            if not self.load_first_sample():
                print("⚠️ Failed to load first sample.")
            
            print("\n📦 Downloading all simulation data...")
            self.download_data()
            
            return True
        except Exception as e:
            print(f"❌ Error in simulation workflow: {e}")
            return False

if __name__ == "__main__":
    # Initialize client with default form data
    client = SimulationClient(
        sequence="SATH",
        samples=10,
        force_field="amber03",
        grid=4,
    )
    
    # Run the simulation workflow
    success = client.run()
    if not success:
        exit(1)