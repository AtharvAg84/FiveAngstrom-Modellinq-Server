import psycopg2
import hashlib
import logging
from typing import Dict, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MoleculeFetcher:
    """Class to fetch molecule data from the database and convert to XYZ format."""
    
    def __init__(self, db_params: Dict[str, str]):
        """
        Initialize with database connection parameters.
        
        Args:
            db_params (dict): Database connection parameters (dbname, user, password, host, port).
        """
        self.db_params = db_params
        self.conn = None
        self.cur = None

    def connect(self) -> None:
        """Establish database connection."""
        try:
            self.conn = psycopg2.connect(**self.db_params)
            self.cur = self.conn.cursor()
            logger.info("Connected to database")
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            raise

    def disconnect(self) -> None:
        """Close database connection."""
        if self.cur:
            self.cur.close()
        if self.conn:
            self.conn.close()
            logger.info("Disconnected from database")

    def compute_composition_hash(self, composition: str) -> Optional[int]:
        """
        Compute a numeric SHA-256 hash for the composition.
        
        Args:
            composition (str): Chemical composition (e.g., 'C9H8O4').
            
        Returns:
            int: 64-bit integer hash, or None if composition is empty.
        """
        if not composition:
            logger.warning("Empty composition provided")
            return None
        hash_bytes = hashlib.sha256(composition.encode('utf-8')).digest()
        hash_int = int.from_bytes(hash_bytes[:8], 'big') & 0x7FFFFFFFFFFFFFFF
        return hash_int

    def query_molecule(self, composition: str) -> Optional[Dict]:
        """
        Query the info table by composition hash.
        
        Args:
            composition (str): Chemical composition.
            
        Returns:
            dict: Molecule data (num_atoms, chemical_formula, chemical_symbol, coordinates, composition),
                  or None if not found.
        """
        logger.info(f"Querying molecule for composition: {composition}")
        hash_int = self.compute_composition_hash(composition)
        if hash_int is None:
            return None

        try:
            self.connect()
            self.cur.execute(
                """
                SELECT num_atoms, chemical_formula, chemical_symbol, coordinates, composition
                FROM info
                WHERE composition_hash = %s
                LIMIT 1
                """,
                (hash_int,)
            )
            result = self.cur.fetchone()
            if not result:
                logger.info("No molecule found for composition hash")
                return None

            num_atoms, chemical_formula, chemical_symbol, coordinates, composition = result
            logger.info(f"Molecule found: {chemical_formula}")
            return {
                "num_atoms": num_atoms,
                "chemical_formula": chemical_formula,
                "chemical_symbol": chemical_symbol,
                "coordinates": coordinates,
                "composition": composition
            }
        except Exception as e:
            logger.error(f"Database query error: {str(e)}")
            return None
        finally:
            self.disconnect()

    def to_xyz(self, data: Dict) -> str:
        """
        Convert molecule data to XYZ format.
        
        Args:
            data (dict): Molecule data with num_atoms, chemical_formula, chemical_symbol, coordinates.
            
        Returns:
            str: XYZ-formatted string.
        """
        if not data or not all(key in data for key in ["num_atoms", "chemical_formula", "chemical_symbol", "coordinates"]):
            logger.warning("Invalid molecule data for XYZ conversion")
            return ""

        num_atoms = data["num_atoms"]
        chemical_formula = data["chemical_formula"]
        chemical_symbol = data["chemical_symbol"]
        coordinates = data["coordinates"]

        if len(chemical_symbol) != num_atoms or len(coordinates) != num_atoms:
            logger.error(f"Mismatch: {len(coordinates)} coordinates, {len(chemical_symbol)} symbols, {num_atoms} atoms")
            return ""

        xyz_lines = [str(num_atoms), chemical_formula]  # Only chemical_formula in comment line
        for symbol, coord in zip(chemical_symbol, coordinates):
            try:
                # Handle list format: [x, y, z]
                if isinstance(coord, list) and len(coord) >= 3:
                    x, y, z = coord[0], coord[1], coord[2]
                # Handle dict format: {"x": x, "y": y, "z": z}
                elif isinstance(coord, dict):
                    x = coord.get("x", 0.0)
                    y = coord.get("y", 0.0)
                    z = coord.get("z", 0.0)
                else:
                    logger.error(f"Invalid coordinate format: {coord}")
                    return ""
                xyz_lines.append(f"{symbol} {x:8.2f} {y:8.2f} {z:8.2f}")
            except Exception as e:
                logger.error(f"Error processing coordinate {coord}: {str(e)}")
                return ""

        return "\n".join(xyz_lines)