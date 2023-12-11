"""
HotWaterDemand_D.py
A. Molar-Cruz @ TUM ENS
"""


class HotWaterDemand_D:
    # --------------------------------------------------------------------------------
    def __init__(self, resolution, day_vector, Tw, daily_DHW):
        """
        Initializes the HotWaterDemand class with the given parameters.

        Args:
            resolution: The resolution of the simulation in minutes.
            day_vector: A list of integers representing the day of the year in the simulation time frame.
            Tw: The supply temperature of water in degrees Celsius.
            daily_DHW: The mean daily hot water consumption in cubic meters.

        """
        # Mean daily consumption
        self.daily_DHW = daily_DHW  # Mean daily hot water consumption [m3]

        # Hot Water values
        self.Tw = Tw  # Supply temperature of water in degC
        self.density = 1000  # Water density in kg/m3
        self.c_p = 4186  # Heat capacity of water in J/(kg degC)

        # Simulation time frame
        self.resolution = resolution  # resolution in min
        self.day_vector = day_vector  # Vector with day of year in simulation time frame

        # Results
        self.dhw_m3 = 0.  # Instantaneous DHW demand in m3
        self.dhw_energy = 0.  # Domestic hot water demand in W

    # --------------------------------------------------------------------------------
    def calculate(self):
        """
        Calculates aggregated hot water demand
        """

        self.dhw_m3 = self.daily_DHW * len(self.day_vector)

        # calculate power consumption in W
        Tw_cold = 6  # degC
        Tw_hot = self.Tw
        dTw = Tw_hot - Tw_cold
        self.dhw_energy = self.dhw_m3 * self.density * self.c_p * dTw / (60 * self.resolution)
