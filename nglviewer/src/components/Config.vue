<template>
    <h3>Configuration</h3>
    <div class="sampling-form">
        <form @submit.prevent="submitForm" ref="configFormRef">
            <!-- Existing Form Fields -->
            <div class="form-field">
                <label for="seq">Sequence</label>
                <div class="aminoSeq">
                    <textarea id="seq" v-model="form.sequence" :class="{ input: true, invalid: invalidSequence }"
                        placeholder="Amino Acid Sequence" rows="1" @input="validateSequence" required
                        :disabled="loading"></textarea>
                    <div v-if="invalidSequence" class="error-message">
                        Invalid Sequence
                    </div>
                </div>
            </div>

            <div class="form-field">
                <label for="samples">Samples</label>
                <input type="number" id="samples" v-model.number="form.samples" class="input"
                    placeholder="No. of Samples" required :disabled="loading" />
            </div>

            <div class="form-field">
                <label for="force">Force Field</label>
                <select id="force" v-model="form.forceField" class="input" :disabled="loading">
                    <option value="amber03.xml">AMBER03</option>
                    <option value="amber10.xml">AMBER10</option>
                    <option value="amber96.xml">AMBER96</option>
                    <option value="amber99sb.xml">AMBER99sb</option>
                    <option value="amberfb15.xml">AMBERFB15</option>
                </select>
            </div>

            <div class="form-field">
                <label for="gridSplit">Grid Split</label>
                <select id="gridSplit" v-model="form.grid" class="input" :disabled="loading">
                    <option value="1">1</option>
                    <option value="2">2</option>
                    <option value="4">4</option>
                    <option value="5">Any</option>
                </select>
            </div>

            <div class="form-field">
                <label for="minimize">Minimize</label>
                <select id="minimize" v-model="form.minimize" class="input" :disabled="loading">
                    <option value="0">None</option>
                    <option value="5">Top 5</option>
                    <option value="10">Top 10</option>
                    <option value="15">Top 15</option>
                    <option value="20">Top 20</option>
                </select>
            </div>

            <div v-if="processId == null" class="button-group">
                <button type="submit" class="generate-button">
                    <i class="fas fa-play"></i>
                    Generate
                </button>
            </div>
        </form>
    </div>
    <div v-if="processId != null" class="progress-container">
        <h4>
            Sampling Progress
            <span class="percentage"> {{ simulationProgress }}%</span>
        </h4>
        <div class="progress-bar">
            <div class="progress-bar-fill" :style="{ width: simulationProgress + '%' }"></div>
        </div>
        <div class="progress-text">
            <div class="status">Status: {{ statusMessage }}</div>
        </div>
    </div>
    <div v-if="simulationProgress == 100" class="button-group">
        <button class="generate-button" @click="downloadData">
            <i class="fas fa-download"></i>
            Download
        </button>
    </div>
</template>

<script>
import axios from "axios";
import JSZip from "jszip";
import { saveAs } from "file-saver";

export default {
    data() {
        return {
            loading: false,
            invalidSequence: false,
            form: {
                sequence: "",
                samples: 100,
                forceField: "amberfb15.xml",
                grid: 4,
                minimize: 0,
            },
            processId: null,
            statusMessage: "Queued",
            simulationProgress: 0,
            statusInterval: null,
            completed: false, // Flag to show download button once completed.
        };
    },
    emits: ["sequence-data", "process-info"],
    methods: {
        validateSequence(e) {
            const allowedChars = "ACDEFGHIKLMNPQRSTVWY";
            let sequence = e.target.value.toUpperCase();
            let isValid = true;
            if (sequence.length === 0) {
                this.invalidSequence = false;
                this.$emit("sequence-data", "");
                return;
            }
            for (let char of sequence) {
                if (!allowedChars.includes(char)) {
                    isValid = false;
                    break;
                }
            }
            this.invalidSequence = !isValid;
            if (isValid) {
                this.form.sequence = sequence;
                this.fetchSequenceData(sequence);
            }
        },
        async submitForm() {
            if (this.invalidSequence) {
                alert("Invalid Sequence");
                return;
            }
            if (this.form.samples < 1 || this.form.samples > 1000) {
                alert("Number of samples must be between 1 and 1000");
                return;
            }
            this.loading = true;
            this.completed = false;

            try {
                const response = await axios.post("/api/sample", {
                    sequence: this.form.sequence,
                    sample: this.form.samples,
                    forcefield: this.form.forceField,
                    grid: this.form.grid,
                    minimize: this.form.minimize,
                });

                console.log(
                    "Config.vue: Response from /api/sample:",
                    response.data,
                );
                this.$emit("process-info", {
                    processId: response.data,
                    totalSamples: this.form.samples,
                    minimize: this.form.minimize,
                });
                console.log("Config.vue: Emitted process-info with:", {
                    processId: response.data,
                    totalSamples: this.form.samples,
                    minimize: this.form.minimize,
                });
                this.processId = response.data;
                this.StatusUpdate();
                console.log("Simulation Data:", this.form);
            } catch (error) {
                console.error("Error starting simulation:", error);
                this.statusMessage = "Error starting simulation";
                this.loading = false;
            }
        },
        fetchSequenceData(sequence) {
            axios
                .get(`/api/seq?${sequence}`)
                .then((response) => {
                    this.$emit("sequence-data", response.data);
                })
                .catch((error) => {
                    console.log("Error fetching sequence data:", sequence);
                    this.invalidSequence = true;
                });
        },
        StatusUpdate() {
            this.statusInterval = setInterval(async () => {
                try {
                    // FIX: Use ?id=... instead of ?...
                    const response = await axios.get(
                        `/api/status?id=${this.processId}`
                    );
                    const newStatus = response.data;
                    this.statusMessage = newStatus;
                    // Update progress based on status
                    if (newStatus.startsWith("Waiting on Queue")) {
                        this.simulationProgress = 0;
                    } else if (newStatus === "Generating the Samples") {
                        this.simulationProgress = 33;
                    } else if (newStatus === "Minimizing the Structures") {
                        this.simulationProgress = 66;
                    } else if (newStatus === "Completed") {
                        this.simulationProgress = 100;
                        this.completed = true; // Setting completed to true when process is completed
                        clearInterval(this.statusInterval);
                        this.loading = false;
                        this.loadFirstSample();
                    }
                } catch (error) {
                    console.error("Error fetching status:", error);
                    this.statusMessage = "Error fetching status";
                    clearInterval(this.statusInterval);
                }
            }, 1000);
        },
        async loadFirstSample() {
            try {
                const sampleUrl = `/api/pdb?Result/${this.processId}/sample/sample_0000.pdb`;
                const response = await axios.get(sampleUrl, {
                    headers: {
                        "Cache-Control": "no-cache",
                        Pragma: "no-cache",
                    },
                });
                if (response.data) {
                    this.$emit("sequence-data", response.data);
                }
            } catch (error) {
                console.error("Error fetching sample data:", error);
            }
        },

        async downloadData() {
            this.loading = true;
            try {
                const zip = new JSZip();
                const processId = this.processId;
                const sampleFolder = zip.folder(processId + "/sample");
                const minimizeFolder = zip.folder(processId + "/minimize");
                // Fetch all the sample pdb data
                for (let i = 0; i < this.form.samples; i++) {
                    const paddedNumber = String(i).padStart(4, "0");
                    const sampleUrl = `/api/pdb?Result/${processId}/sample/sample_${paddedNumber}.pdb`;
                    const sampleOutUrl = `/api/pdb?Result/${processId}/sample/sampled.out`;

                    const [sampleResponse, sampleOutResponse] =
                        await Promise.all([
                            axios.get(sampleUrl),
                            axios.get(sampleOutUrl),
                        ]);
                    if (sampleResponse.data) {
                        sampleFolder.file(
                            `sample_${paddedNumber}.pdb`,
                            sampleResponse.data,
                        );
                    }
                    if (sampleOutResponse.data) {
                        sampleFolder.file(
                            "sampled.out",
                            sampleOutResponse.data,
                        );
                    }
                }
                // Fetch the minimized data
                if (this.form.minimize != 0) {
                    for (let i = 0; i < this.form.minimize; i++) {
                        const paddedNumber = String(i).padStart(4, "0");
                        const minimizeUrl = `/api/pdb?Result/${processId}/minimize/minimized_${paddedNumber}.pdb`;
                        const minimizedOutUrl = `/api/pdb?Result/${processId}/minimize/minimized.out`;

                        const [minimizeResponse, minimizedOutResponse] =
                            await Promise.all([
                                axios.get(minimizeUrl),
                                axios.get(minimizedOutUrl),
                            ]);
                        if (minimizeResponse.data) {
                            minimizeFolder.file(
                                `minimized_${paddedNumber}.pdb`,
                                minimizeResponse.data,
                            );
                        }
                        if (minimizedOutResponse.data) {
                            minimizeFolder.file(
                                "minimized.out",
                                minimizedOutResponse.data,
                            );
                        }
                    }
                }

                const zipBlob = await zip.generateAsync({ type: "blob" });
                saveAs(zipBlob, `${processId}.zip`); // Changed the file name here
            } catch (error) {
                console.error("Error downloading data:", error);
            } finally {
                this.loading = false;
            }
        },
    },
};
</script>

<style scoped>
.sampling-form {
    width: 100%;
    box-sizing: border-box;
    padding: 1rem;
}

.form-field {
    display: flex;
    margin-bottom: 1rem;
    width: 100%;
    box-sizing: border-box;
}

.form-field label {
    margin-right: 1rem;
    font-weight: 900;
    color: var(--text-color);
    font-size: 0.9rem;
    white-space: nowrap;
    width: 100px;
    align-self: center;
}

.input {
    width: 100%;
    box-sizing: border-box;
    padding: 0.6rem;
    background-color: var(--background-color);
    border: 1px solid var(--border-color);
    border-radius: 4px;
    font-family: "Helvetica Neue", Arial, "Segoe UI", Roboto, "Open Sans";
    font-size: 0.9rem;
    color: var(--text-color);
    transition: all 0.2s ease;
}

.input:focus {
    outline: none;
    border-color: var(--primary-color, #4caf50);
    box-shadow: 0 0 0 1px var(--secondary-color, #4caf50);
}

.input.invalid {
    border-color: var(--error-color, #f44336);
}

textarea.input {
    resize: vertical;
    min-height: 60px;
}

.input-checkbox {
    width: auto;
    margin-left: 10px;
    align-self: center;
    box-sizing: border-box;
    padding: 0.6rem;
    background-color: var(--background-color);
    border: 1px solid var(--border-color);
    border-radius: 4px;
    transition: all 0.2s ease;
}

.generate-button {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 0.5rem;
    width: 100%;
    padding: 0.8rem 1.5rem;
    background: linear-gradient(135deg,
            var(--primary-color),
            var(--secondary-color));
    color: white;
    border: none;
    border-radius: 25px;
    font-weight: 500;
    font-size: 1rem;
    cursor: pointer;
    transition: all 0.3s ease;
    box-shadow: 0 2px 5px rgba(0, 0, 0, 0.2);
}

.generate-button:hover:not(:disabled) {
    transform: translateY(-1px);
    box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
}

.generate-button:disabled {
    background: var(--disabled-color);
    cursor: not-allowed;
    transform: none;
    box-shadow: none;
}

.generate-button i {
    font-size: 0.9em;
}

.button-group {
    display: flex;
    gap: 0.8rem;
    margin-top: 1.5rem;
    width: 100%;
    box-sizing: border-box;
}

.btn {
    flex: 1;
    padding: 0.7rem 1rem;
    border: none;
    border-radius: 4px;
    font-weight: 500;
    font-size: 0.9rem;
    cursor: pointer;
    transition: all 0.2s ease;
    white-space: nowrap;
}

.submit {
    background-color: var(--accent-color, #4caf50);
    color: white;
}

.submit:hover:not(:disabled) {
    opacity: 0.9;
}

.submit:disabled {
    background-color: var(--disabled-color, #cccccc);
    cursor: not-allowed;
}

.cancel {
    background-color: var(--error-color, #f44336);
    color: white;
}

.cancel:hover {
    opacity: 0.9;
}

.aminoSeq {
    display: flex;
    flex-direction: column;
}

.error-message {
    color: var(--error-color, #f44336);
    font-size: 0.8rem;
    display: block;
    text-align: center;
}

.progress-container {
    margin-top: 1rem;
    padding: 1rem;
    border-top: 1px solid var(--border-color);
    display: flex;
    flex-direction: column;
    align-items: flex-start;
}

.progress-container h4 {
    display: flex;
    align-items: center;
    margin-bottom: 0.3rem;
}

.progress-container .percentage {
    margin-left: 0.5rem;
    font-weight: bold;
    font-size: 1rem;
    color: var(--text-color);
}

.progress-bar {
    width: 100%;
    height: 15px;
    background-color: #e0e0e0;
    border-radius: 7px;
    overflow: hidden;
    margin-top: 0.5rem;
    margin-bottom: 0.5rem;
}

.progress-bar-fill {
    height: 100%;
    background: linear-gradient(to right,
            var(--primary-color),
            var(--secondary-color));
    border-radius: 7px;
    transition: width 0.3s ease;
}

.progress-text {
    font-size: 0.9rem;
    color: var(--text-color-secondary);
    text-align: center;
    font-style: italic;
}

@media (max-width: 400px) {
    .sampling-form {
        padding: 0.8rem;
    }

    .button-group {
        flex-direction: column;
    }

    .btn {
        width: 100%;
    }
}
</style>
