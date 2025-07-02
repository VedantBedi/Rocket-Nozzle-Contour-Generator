# Rocket Nozzle Generator

This Python project generates rocket nozzle contours using the Method of Characteristics (MOC). It calculates the throat diameter from mass flow rate and propulsion inputs, creates a nozzle profile, and exports contour data to a CSV file. This also has a dusion script so that we can get the model straight in Fusion360.

## Features
- Calculates throat diameter
- Generates nozzle profile using MOC
- Exports contour points to CSV
- Plots Mach number distribution

## Installation
Clone the repository:
```
git clone https://github.com/EakamjitSingh/rocket-nozzle-generator
cd rocket-nozzle-generator
```
Install dependencies:
```
pip install -r requirements.txt
```

## Usage
Run the script:
```
python src/rocket_contour_generator.py
```
Provide inputs when prompted:
- Chamber Pressure (Pa)
- Chamber Temperature (K)
- Mass Flow Rate (kg/s)
- Exit Mach Number

The program generates the nozzle contour, plots Mach distribution, and exports data to `nozzle_contour_points.csv`.

## Output
- Nozzle contour plot
- CSV with (x, y, z) coordinates

## Project Structure
```
rocket-nozzle-generator/
│
├── src/                      # Source code
│   └── rocket_contour_generator.py
|   └── rocket_contour_fusion_cript.py
│
|
│
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```

## Contributing
Fork, create a branch, and submit a pull request.

## License
MIT License.
