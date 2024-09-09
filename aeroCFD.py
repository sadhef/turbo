# Class defining aero CFD generation framework.

# Import statements.
#import trieste as tst
#import tensorflow as tf
#import pandas as pd
import json
import os


# WHAT:   Class of controls for iteratiely optimised CFD dataset generation. 
# WHY:    Avoid generating low-value-add CFD data (e.g. via uniform parameter distrubution). 
# OUTPUT: N/A.
class AeroCFD:

    def __init__(self):
        pass

    # WHAT:   Generate CFD dataset, guided iteratively by a Deep GPR model find high-value data.
    # WHY:    Avoid generating low value data, ranked by eta_uniformJ performance metric.
    # OUTPUT: CFD performance dataset, ranked by eta_uniformJ performance metric.
    def generate_iterative_dataset(self, cfds_per_iteration: int, limits: json):
        # Set observer for CFD performance data point evaluations of given design (d),
        # based on a asychnronous evaluator (asynchronous due to CFD generation times and queues).

        pass

    # WHAT:   Request HPC to generate specified CFDs to evaluate given designs.
    # WHY:    To use CFDs as objective function's evaluator, for low quasi-parallel compute. 
    # OUTPUT: 3D Array - Float : CFD performance results for given designs.  
    def objective(self, designs_for_evaluation):

        dummy_data = {
            "name": "Eryk Krusinski",
            "email": "erykkrusinski@example.com",
            "courses": [
                {
                    "course_name": "Machine Learning",
                    "grade": "A"
                },
                {
                    "course_name": "Data Structures",
                    "grade": "A"
                }
            ]

        }

        # Create path of current directory and filename.
        filepath = os.path.join(os.getcwd(), "input_design_data.json")

        # Write the JSON data to the file.
        with open(filepath, 'w') as json_file:
            json.dump(dummy_data, json_file, indent=4)

        pass
