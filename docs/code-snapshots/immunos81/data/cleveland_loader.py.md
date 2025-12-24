---
source: /Users/byron/projects/immunos81/data/cleveland_loader.py
relative: immunos81/data/cleveland_loader.py
generated_at: 2025-12-23 10:28
---

```python
"""
Cleveland Heart Disease Dataset Loader for Immunos-81

The Cleveland dataset contains 303 patient records with 14 attributes used to
diagnose coronary artery disease (CAD).
"""

from typing import List, Tuple
import numpy as np
from ..core.amino_acid_library import AminoAcidLibrary
from ..core.antigen import Antigen


class ClevelandDataLoader:
    """
    Loader for the Cleveland Heart Disease dataset.

    The dataset contains 303 instances with 14 attributes:
    1. age: age in years (numeric)
    2. sex: sex (1 = male, 0 = female) (nominal)
    3. cp: chest pain type (nominal, 1-4)
    4. trestbps: resting blood pressure (numeric)
    5. chol: serum cholesterol in mg/dl (numeric)
    6. fbs: fasting blood sugar > 120 mg/dl (nominal, 0/1)
    7. restecg: resting electrocardiographic results (nominal, 0-2)
    8. thalach: maximum heart rate achieved (numeric)
    9. exang: exercise induced angina (nominal, 0/1)
    10. oldpeak: ST depression induced by exercise (numeric)
    11. slope: slope of peak exercise ST segment (nominal, 1-3)
    12. ca: number of major vessels colored by fluoroscopy (numeric, 0-3)
    13. thal: thalassemia (nominal, 3/6/7)
    14. num: diagnosis (0 = no disease, 1-4 = disease severity)

    For Immunos-81, we binarize the diagnosis: CAD- (no disease) vs CAD+ (disease present)
    """

    # Cleveland dataset URL from UCI repository
    DATA_URL = "https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data"

    def __init__(self):
        """Initialize the dataset loader"""
        self.library = AminoAcidLibrary()
        self._setup_library()

    def _setup_library(self):
        """Set up the Amino Acid Library with all dataset variables"""
        # Add all 13 input variables (excluding the target)
        self.library.add_variable("age", "numeric")
        self.library.add_variable("sex", "nominal")
        self.library.add_variable("cp", "nominal")
        self.library.add_variable("trestbps", "numeric")
        self.library.add_variable("chol", "numeric")
        self.library.add_variable("fbs", "nominal")
        self.library.add_variable("restecg", "nominal")
        self.library.add_variable("thalach", "numeric")
        self.library.add_variable("exang", "nominal")
        self.library.add_variable("oldpeak", "numeric")
        self.library.add_variable("slope", "nominal")
        self.library.add_variable("ca", "numeric")
        self.library.add_variable("thal", "nominal")

    def load_data(self) -> Tuple[List[Antigen], AminoAcidLibrary]:
        """
        Load the Cleveland dataset.

        Returns:
            Tuple of (list of antigens, amino acid library)

        Note:
            This attempts to download from UCI repository. If that fails,
            you should provide the data file locally.
        """
        try:
            import urllib.request
            print(f"Downloading Cleveland dataset from UCI repository...")
            response = urllib.request.urlopen(self.DATA_URL)
            data_str = response.read().decode('utf-8')
            lines = data_str.strip().split('\n')
            print(f"Downloaded {len(lines)} records")
        except Exception as e:
            print(f"Could not download dataset: {e}")
            print("Generating synthetic demo data instead...")
            return self._generate_demo_data()

        antigens = []
        for i, line in enumerate(lines):
            if not line.strip():
                continue

            try:
                antigen = self._parse_line(line, i)
                if antigen is not None:
                    antigens.append(antigen)
            except Exception as e:
                print(f"Warning: Could not parse line {i}: {e}")
                continue

        print(f"Successfully loaded {len(antigens)} antigens")
        return antigens, self.library

    def _parse_line(self, line: str, index: int) -> Antigen:
        """Parse a single line from the Cleveland dataset"""
        values = line.strip().split(',')

        if len(values) != 14:
            return None

        # Handle missing values (marked as '?')
        processed_values = []
        for val in values[:-1]:  # Exclude target variable
            if val == '?':
                processed_values.append(None)
            else:
                try:
                    # Try to convert to float
                    processed_values.append(float(val))
                except ValueError:
                    processed_values.append(val)

        # Get target value and convert to binary class
        try:
            target = int(float(values[-1]))
            class_label = "CAD-" if target == 0 else "CAD+"
        except:
            return None

        # Create antigen
        variable_names = ["age", "sex", "cp", "trestbps", "chol", "fbs",
                         "restecg", "thalach", "exang", "oldpeak",
                         "slope", "ca", "thal"]

        data_dict = {name: val for name, val in zip(variable_names, processed_values)}

        antigen = Antigen.from_dict(
            data_dict,
            self.library,
            class_label=class_label,
            identifier=f"cleveland_{index}"
        )

        return antigen

    def _generate_demo_data(self) -> Tuple[List[Antigen], AminoAcidLibrary]:
        """
        Generate synthetic demo data for testing.

        This creates realistic-looking data with similar distributions
        to the Cleveland dataset.
        """
        np.random.seed(42)
        antigens = []

        n_samples = 303

        for i in range(n_samples):
            # Generate synthetic patient data
            age = np.random.normal(54, 9)  # Mean age ~54
            sex = np.random.choice([0, 1])
            cp = np.random.choice([1, 2, 3, 4])
            trestbps = np.random.normal(131, 17)
            chol = np.random.normal(246, 51)
            fbs = np.random.choice([0, 1], p=[0.85, 0.15])
            restecg = np.random.choice([0, 1, 2], p=[0.5, 0.45, 0.05])
            thalach = np.random.normal(150, 23)
            exang = np.random.choice([0, 1], p=[0.67, 0.33])
            oldpeak = np.random.exponential(1.0)
            slope = np.random.choice([1, 2, 3], p=[0.46, 0.41, 0.13])
            ca = np.random.choice([0, 1, 2, 3], p=[0.59, 0.21, 0.12, 0.08])
            thal = np.random.choice([3, 6, 7], p=[0.03, 0.51, 0.46])

            # Class label (roughly 54% disease, 46% no disease)
            has_disease = np.random.random() < 0.54
            class_label = "CAD+" if has_disease else "CAD-"

            data_dict = {
                "age": age,
                "sex": sex,
                "cp": cp,
                "trestbps": trestbps,
                "chol": chol,
                "fbs": fbs,
                "restecg": restecg,
                "thalach": thalach,
                "exang": exang,
                "oldpeak": oldpeak,
                "slope": slope,
                "ca": ca,
                "thal": thal
            }

            antigen = Antigen.from_dict(
                data_dict,
                self.library,
                class_label=class_label,
                identifier=f"demo_{i}"
            )
            antigens.append(antigen)

        print(f"Generated {len(antigens)} synthetic demo records")
        return antigens, self.library

```
