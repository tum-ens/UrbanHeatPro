"""
SpaceHeatingDemand.py
A. Molar-Cruz @ TUM ENS
"""

import copy
from random import randint

import numpy as np


# import ipdb


class SpaceHeatingDemand:
    # --------------------------------------------------------------------------------
    def __init__(self, dt_vector, resolution, heated_area,
                 Tamb, I, Tb0, dT_per_hour, eta,
                 thermal_intertia, U, V, C, Tset, dTset,
                 activity_vector, occupancy_vector, sh_prob,
                 _solar_gains, _internal_gains,
                 _night_set_back, schedule_nsb, T_nsb, power_reduction,
                 window_areas, coords, debug):
        """
        Initializes an instance of the SpaceHeatingDemand class.

        Args:
            dt_vector: List of time steps as datetime objects.
            resolution: Resolution in minutes.
            heated_area: Heated area in square meters.
            Tamb: Ambient temperature vector in degrees Celsius.
            I: Solar radiation vector in W/m2.
            Tb0: Initial building temperature in degrees Celsius.
            dT_per_hour: Maximum change in temperature allowed per hour in degrees Celsius.
            eta: Heating process efficiency.
            thermal_intertia: Thermal inertia of the heating system.
            U: Building transmission losses in W/K.
            V : Building ventilation losses in W/K.
            C: Equivalent building thermal mass in J/K.
            Tset: Set temperature or target temperature in degrees Celsius.
            dTset: Delta temperature for Tset_min and Tset_max.
            activity_vector: Building activity vector (0, 1).
            occupancy_vector: Number of occupants in the building in each time step.
            sh_prob: Probability vector of using space heating.
            _solar_gains: Solar gains in W/m2.
            _internal_gains: Internal gains in W/m2.
            _night_set_back: Share of buildings with night set-back.
            schedule_nsb: Start and end of night set-back in hours.
            T_nsb: Night set-back temperature in degrees Celsius.
            power_reduction: Percentage of power reduced (as decimal).
            window_areas: Window area oriented to [E, S, W, N] in square meters.
            coords: (latitude, longitude) of the building centroid.
            debug: Debug flag.
        """

        # Building thermal properties
        self.U = U  # Building transmission losses in W/K
        self.V = V  # Building ventilation losses in W/K
        self.C = C  # Equivalent building thermal mass in J/K
        self.heated_area = heated_area  # Heated area in m2

        # Heat gains
        self._solar_gains = _solar_gains
        self._internal_gains = _internal_gains
        self.window_areas = window_areas  # Window area oriented to [E, S, W, N] in m2
        self.coords = coords  # (lat, lon) of building centroid

        # User behavior
        self.sh_prob = sh_prob  # Probability vector of using space heating
        self.activity_vector = activity_vector  # Building activity vector (0, 1)
        self.occupancy_vector = occupancy_vector  # Number of occupants in building in each time step
        self.Tset = Tset  # Set temperature or target temperature in degC
        self.Tset_original = copy.deepcopy(self.Tset)  # Copy of original set temperature in degC
        self.dTset = dTset  # Delta temperature (for Tset_min, Tset_max)
        self.Tb0 = Tb0  # Initial building temperature (Tb @ t=0) in degC

        # External factors
        self.Tamb = Tamb  # Ambient temperature vector in degC
        self.I = I  # Solar radiation vector in W/m2 [I_Gh, I_Dh, I_ex, hs]
        self.eta = eta  # Heating process efficiency
        self.thermal_intertia = thermal_intertia  # Thermal inertia of the heating system

        # Maximum power
        self.dT_per_hour = dT_per_hour  # Maximum dT allowed in building per hour [degC]
        self.sh_powermax = 0.  # Maximum power per timestep

        # DSM
        self._night_set_back = _night_set_back  # Share of buildings with nsb
        self.schedule_nsb = schedule_nsb  # [start, end] of nsb in h
        self.T_nsb = T_nsb  # Night set-back temperature in degC
        self.power_reduction = power_reduction  # Percentage of power reduced (as decimal)

        # Simulation time steps
        self.dt_vector = dt_vector  # Vector of time steps as datetime objects
        self.nts = len(self.dt_vector)  # Number of time steps
        self.resolution = resolution  # Resolution in min

        # Time series
        self.sh_power = np.zeros([self.nts])  # Space heating demand W
        self.internal_gains = np.zeros([self.nts])  # Internal gains in W
        self.solar_gains = np.zeros([self.nts])  # Solar gains in W
        self.Tb = np.zeros([self.nts])  # Building temperature in degC

        #
        self.debug = debug

    # --------------------------------------------------------------------------------
    def calculate(self):
        """
        Calculates the time series of space heating demand for a single building
        as the numerical solution of a first order building thermal model (1R1C).
        Transmission and ventilation losses through infiltration are included.

        returns:
            self.Tb
            self.sh_power
            self.internal_gains
            self.solar_gains
        """

        # Initial conditions
        self.Tb[0] = self.Tb0
        # print(0,self.sh_power[0], self.Tb[0])

        # Solar gains
        if self._solar_gains:
            # initialize parameters
            ### Reduction factors from TABULA
            F_sh = 0.7  # reduction factor due to external shading
            F_F = 0.25  # fraction of the window area that is frame
            F_W = 1.0  # reduction factor to consider radiation non-perpendicular to the glazing 		>>> check if required
            g_gl = 0.85  # total solar energy transmittance for radiation perpendicular to the glazing
            RED_FACTORS = [F_sh, F_F, F_W, g_gl]

            ### recalculate window orientation
            # [E, S, W, N] -> [-90, 0, 90, -180]
            deviation_from_south = 0
            ORIENTATION = [-90 + deviation_from_south, 0 + deviation_from_south, 90 + deviation_from_south,
                           -180 + deviation_from_south]

        # Heat demand calculation for each time step
        for iii in range(1, self.nts):

            # if self.debug == 1:
            # Progress bar
            # sys.stdout.write('\r')
            # i = float(iii + 1) / self.nts
            # sys.stdout.write("        Space heating: [%-20s] %d%%" % ('|' * int(20 * i), 100 * i))
            # sys.stdout.flush()

            # clean Tamb data
            # if value is inf take value from previous time step
            if np.isinf(self.Tamb[iii]):
                self.Tamb[iii] = self.Tamb[iii - 1]

            # calculate Tset
            self.calculate_Tset(iii)

            # calculate active flag (is heating system on/off?)
            # [temperature_flag, activity_flag, active_flag]
            active_flag = self.calculate_flags(iii)

            # calculate heat gains

            ### calculate heat gains only when heating is required
            ### ONLY FOR COMPARISON (TABULA VS UrbanHeatPro)
            if active_flag[0]:
                if self._internal_gains:
                    self.calculate_internal_gains(iii)
                if self._solar_gains:
                    self.calculate_solar_gains(iii, RED_FACTORS, ORIENTATION)
            # if self._internal_gains:
            #	self.calculate_internal_gains(iii)
            # if self._solar_gains:
            #	self.calculate_solar_gains(iii, RED_FACTORS, ORIENTATION)

            # get time delta in seconds
            dt = self.dt_vector[iii] - self.dt_vector[iii - 1]
            dt_seconds = dt.seconds

            # calculate maximum power based on max delta temperature in building
            dT_in_resolution = self.dT_per_hour / 60. * self.resolution
            self.sh_powermax = (self.C / self.eta) * ((dT_in_resolution / dt_seconds) -
                                                      ((self.U + self.V) / self.C) *
                                                      (self.Tamb[iii - 1] - self.Tb[iii - 1]) -
                                                      (1. / self.C) * (self.solar_gains[iii - 1] + self.internal_gains[
                        iii - 1]))

            # no cooling
            if self.sh_powermax < 0:
                self.sh_powermax = 0

            # First order thermal model
            ### Power
            self.sh_power[iii] = (self.C / self.eta) * (((self.Tset - self.Tb[iii - 1]) / dt_seconds) -
                                                        ((self.U + self.V) / self.C) *
                                                        (self.Tamb[iii - 1] - self.Tb[iii - 1]) -
                                                        (1. / self.C) * (
                                                                    self.solar_gains[iii - 1] + self.internal_gains[
                                                                iii - 1])) * \
                                 active_flag[2] * (1 - self.power_reduction)
            # no cooling
            if self.sh_power[iii] < 0:
                self.sh_power[iii] = 0

            # limit input power
            if iii > 0:
                self.sh_power[iii] = (1. - self.thermal_intertia) * min(self.sh_power[iii], self.sh_powermax) + \
                                     self.sh_power[iii - 1] * self.thermal_intertia
            else:
                self.sh_power[iii] = min(self.sh_power[iii], self.sh_powermax)

            ### Building temperature
            self.Tb[iii] = self.Tb[iii - 1] + dt_seconds * \
                           (1. / self.C) * ((self.U + self.V) * (self.Tamb[iii - 1] - self.Tb[iii - 1]) +
                                            self.eta * self.sh_power[iii - 1] +
                                            self.solar_gains[iii - 1] + self.internal_gains[iii - 1])

    #
    def calculate_Tset(self, iii):
        """
        Returns Tset to original value or recalculates it depending on night set-back
        """

        # set Tset to input original value
        self.Tset = self.Tset_original

        # Night set-back
        # From 23:00 to 6:00 temperature is lowered to 18degC
        # check if building has night set-back
        rand_num = np.random.uniform(0, 1, 1)
        if rand_num < self._night_set_back:

            # check night-set-back hours
            if self.dt_vector[iii].hour >= self.schedule_nsb[0] or self.dt_vector[iii].hour < self.schedule_nsb[1]:
                # set Tset at Tnight_set_back
                self.Tset = self.T_nsb

    #
    def calculate_flags(self, iii):
        """
        Calculates if heating system is active based on building temperature and building occupancy
        """

        temperature_flag = False
        activity_flag = False
        active_flag = False

        ### Probability of using space heating
        rand_num = np.random.uniform(0, 1, 1)
        if rand_num <= self.sh_prob[iii]:

            ### Confort zone
            #   if heating is off, it remains off until Tset_min
            #   if heating is on, it remains on til Tset_max
            # calculate Tset [min, max]
            Tset_min = self.Tset - self.dTset
            Tset_max = self.Tset + self.dTset
            if self.sh_power[iii - 1] == 0:  # off
                # Building temperature below Tset_min
                temperature_flag = self.Tb[iii - 1] < Tset_min
            else:  # on
                # Building temperature below Tset_max
                temperature_flag = self.Tb[iii - 1] < Tset_max

            ### No building activity at time step
            activity_flag = self.activity_vector[iii]

            ### is heat generation active?
            active_flag = float(temperature_flag * activity_flag)

        return [temperature_flag, activity_flag, active_flag]

    #
    def calculate_internal_gains(self, iii):
        """
        Calculates heat gain in time step due to the activeness of the occupants:
            - 80 W/occupant during the night (23:00 to 6:00)
            - Random between 100 - 125 W/occupant for the rest of the day
        From Validation of RC Building Models for Application in Energy and DSM (Kuniyoshi, 2017)
        EESC Kramer
        [VDI 2078]

        returns:
            float:		self.internal_gains[iii]: Heat gain in W
        """

        # Internal heat gain at night
        if self.dt_vector[iii].hour <= 6:
            hg_per_occupant = 80  # in W
        else:
            hg_per_occupant = randint(100, 125)  # in W

        # Calculate total heat gain (internal)
        ### Calculate total heat gain (internal) adjusted to TABULA values of 16 kWh/(m2 a)
        ### ONLY FOR COMPARISON (TABULA vs UrbanHeatPro)
        # self.internal_gains[iii] = 6.25 * self.heated_area

        ### Heat gains from Kuniyoshi 2017
        self.internal_gains[iii] = self.occupancy_vector[iii] * hg_per_occupant * 2

    #
    def calculate_solar_gains(self, iii, RED_FACTORS, ORIENTATION):
        """
        Calculates solar gains based on the windows size and orientation.
        Method adapted from TABULA

        returns:
            float:		self.solar_gains[iii]: Heat gain in W
        """

        # calculate global incident solar radiation on tilted surface in W
        day_of_year = self.dt_vector[iii].timetuple().tm_yday
        hour = self.dt_vector[iii].hour
        I_Gh = self.I[iii][0]
        I_Dh = self.I[iii][1]
        I_ex = self.I[iii][2]
        hs = self.I[iii][3]
        lat = self.coords[0]
        lon = self.coords[1]
        slope = 90.  # vertical windows

        for window in range(4):  # windows are distributed in four orientations only
            # [E, S, W, N]
            orientation = ORIENTATION[window]

            # calculate incident solar radiation in W/m2
            I_Gt = self.calculate_incident_solar_irradiation(day_of_year, hour, I_Gh, I_Dh, I_ex, hs, lat, lon, slope,
                                                             orientation)

            # calculate window area in specific orientation
            area = self.window_areas[window]

            # calculate solar heat gain
            self.solar_gains[iii] += RED_FACTORS[0] * (1 - RED_FACTORS[1]) * RED_FACTORS[2] * RED_FACTORS[
                3] * area * I_Gt

    #
    def calculate_incident_solar_irradiation(self, day_of_year, hour, I_Gh, I_Dh, I_ex, hs, lat, lon, slope,
                                             orientation):
        """
        Calculates the global incident solar irradiation on tilted surface in W/m2.
        Based on HDKR radiation model (anisotropic model) from
        High-resolution spatio-temporal matching of residential electricity consumption and
        PV power production at the urban scale (Molar-Cruz, 2015)

        Args:
            I_Gh		(float):		Global horizonal radiation in W/m2
            I_Dh		(float):		Diffuse horizontal radiation in W/m2
            I_ex		(float):		Extraterrestrial solar radiation in W/m2
            hs			(float):		Sun elevation angle in deg
            lat			(float):		Latitude in degrees
            lon			(float):		Longitude in degrees
            slope		(int):		Inclination angle of window. Vertical = 90 deg
            orientation (int):		Window orientation

        returns:
            float:		I_Gt: Incident global solar radiation on tilted surface
        """

        # check if there is sunlight
        if I_Gh > 0:

            # Declination angle  [deg]
            declination_angle = 23.45 * np.sin((360. / 365.) * (284 + day_of_year) * np.pi / 180.)

            # Time zone east of GMT [h]
            Z_c = 1.

            # Equation of time [h]
            B = (day_of_year - 1) * 360. / 365.
            E = 3.82 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B) - 0.014615 * np.cos(
                2 * B) - 0.04089 * np.sin(2 * B))

            # Solar time [g]
            t_s = hour + (lon / 15.) - Z_c + E

            # Hour angle [deg]
            hour_angle = (t_s - 12) * 15

            # Angle of incidence [deg]
            incidence_angle = np.arccos(
                np.sin(declination_angle * np.pi / 180.) * np.sin(lat * np.pi / 180.) * np.cos(slope * np.pi / 180.) - \
                np.sin(declination_angle * np.pi / 180.) * np.cos(lat * np.pi / 180.) * np.sin(
                    slope * np.pi / 180.) * np.cos(orientation * np.pi / 180.) + \
                np.cos(declination_angle * np.pi / 180.) * np.cos(lat * np.pi / 180.) * np.cos(
                    slope * np.pi / 180.) * np.cos(hour_angle * np.pi / 180.) + \
                np.cos(declination_angle * np.pi / 180.) * np.sin(lat * np.pi / 180.) * np.sin(
                    slope * np.pi / 180.) * np.cos(orientation * np.pi / 180.) * np.cos(hour_angle * np.pi / 180.) + \
                np.cos(declination_angle * np.pi / 180.) * np.sin(slope * np.pi / 180.) * np.sin(
                    orientation * np.pi / 180.) * np.sin(hour_angle * np.pi / 180.)) * (180. / np.pi)

            # zenith_angle
            zenith_angle = 90. - hs

            # Geometric factor R_B
            # Ratio of beam radiation on the tilted surface to beam radiation on the horizontal surface
            if zenith_angle < 89.:
                R_B = np.cos(incidence_angle * np.pi / 180.) / np.cos(zenith_angle * np.pi / 180.)
                if R_B < 0:
                    R_B = 0
            else:
                R_B = 0

            # Beam radiation on tilted surface, I_Bt [W/m2]
            # If the incidence angles is greater than 90deg, it means the sun is
            # behind the surface and therefore there is no direct or beam radiation
            if incidence_angle > 90.:
                I_Bt = 0
            else:
                I_Bt = abs((I_Gh - I_Dh)) * R_B

            # Anisotropy index - Function of the transmittance of the atmosphere for beam radiation.
            # Portion of the horizontal diffuse to be treated as forward scattered.
            if I_ex > I_Dh:
                A_i = I_Dh / I_ex
            else:
                A_i = 1

            # Modulating function
            # Function to become an all sky-model (Klucher model) - Used to account for horizontal brightening
            # and related to cloudiness
            F = np.sqrt(abs((I_Gh - I_Dh)) / I_Gh)

            # Diffuse radiation (tilted), I_Dt [W/m2]
            if incidence_angle > 90:
                I_Dt = I_Dh * (1 - A_i) * (0.5 * (1 + np.cos(slope * np.pi / 180.))) * (
                            1 + F * ((np.sin(slope / 2 * np.pi / 180.)) ** 3.))
            else:
                I_Dt = I_Dh * (1 - A_i) * (0.5 * (1 + np.cos(slope * np.pi / 180.))) * (
                            1 + F * ((np.sin(slope / 2 * np.pi / 180.)) ** 3.)) + A_i * R_B

            # Reflected radiation (tilted), I_Rt [W/m2]
            albedo = 0.2
            # Nominal value for green vegetation and some soil types
            # From: http://rredc.nrel.gov/solar/pubs/bluebook/appendix.html
            I_Rt = 0.5 * albedo * I_Gh * (1 - np.cos(slope * np.pi / 180.))

            # Incident global solar radiation (tilted), G_Gt [W/m2]
            I_Gt = I_Bt + I_Dt + I_Rt

        # no sunlight
        else:
            I_Gt = 0

        if I_Gt < 0:
            print('I_Bt + I_Dt + I_Rt')
            print(I_Bt, I_Dt, I_Rt)
            print(A_i)
            # this is a bug: where does this function come from?
            raw_input()

        return I_Gt
