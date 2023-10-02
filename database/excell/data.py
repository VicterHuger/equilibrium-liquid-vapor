TEMPERATURE = "temperature"
PRESSURE = "pressure"

Fluid1 = [
    {TEMPERATURE: 231.8, PRESSURE: 87.89439},
    {TEMPERATURE: 239.9, PRESSURE: 88.43218},
    {TEMPERATURE: 249.5, PRESSURE: 85.23991},
    {TEMPERATURE: 257.2, PRESSURE: 78.46925},
    {TEMPERATURE: 266.1, PRESSURE: 31.39873},
    {TEMPERATURE: 267.1, PRESSURE: 57.59882},
]

Fluid2 = [
    {TEMPERATURE: 235.99, PRESSURE: 96.2508},
    {TEMPERATURE: 247.42, PRESSURE: 103.3937},
    {TEMPERATURE: 258.88, PRESSURE: 106.5791},
    {TEMPERATURE: 270.14, PRESSURE: 106.324},
    {TEMPERATURE: 289.5, PRESSURE: 96.80238},
    {TEMPERATURE: 294.03, PRESSURE: 12.99455},
    {TEMPERATURE: 294.97, PRESSURE: 88.55625},
    {TEMPERATURE: 298.64, PRESSURE: 80.13776},
    {TEMPERATURE: 303.71, PRESSURE: 64.26948},
    {TEMPERATURE: 304.52, PRESSURE: 31.69244},
]


def create_fluid_data(name: str):
    """Function to create fluid data to database insertion

    Args:
        name (str): name of the fluid that want to create data, it is necessary to have a list with this name

    Returns:
        _type_: return the data of insertion in database 
    """
    return {
        "name": name,
        "temp_pressure_coordinates": [{TEMPERATURE: item[TEMPERATURE], PRESSURE:item[PRESSURE]} for item in globals()[name]]
    }
