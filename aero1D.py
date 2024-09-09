"""Class defining 1-D Aero model"""

# import statements
import numpy as np
import matplotlib.pyplot as plt
import json
import copy

class Aero1D:
    """Defines 1-D Aero model
    
    Parameters
    ----------
    aero_input : string or dict, optional
        Dictionary or path to the JSON object specifying aero1D input parameters
        
        name : str, optional
            Name of the aerodynamic model
        
        turbomachinary: dictionary
            Definition of turbomachinery parameters
        
            tip_dia : float
                Tip diameter of the rotor in m
            
            htr : float, optional
                Hub to tip ratio
            
            hub_dia : float
                Hub diameter of the rotor in m
            
            sigma : float, optional
                Area ratio
        
        design-mean-line: dictionary
            Definition of design mean line parameters
            
            phi : float, optional
                Flow coefficient
            
            eta_turbo : float, optional
                Fan efficiency. The default is 0.8
    
    """
    def __init__(self, aero_input):
        
        # File
        if isinstance(aero_input, str):
            with open(aero_input) as input_json_handle:
                self._input_dict = json.load(input_json_handle)
        
        # Dictionary
        elif isinstance(aero_input, dict):
            self._input_dict = copy.deepcopy(aero_input)
        
        else:
            raise IOError("Input to the class initialisation must be a file or python dictionary, not type {0}.".format(type(aero_input)))
        
        self.name = self._input_dict.get("name", "Unnamed")
        
        _fan_geom_dict = self._input_dict.get("turbomachinery",{}).get("fan-geometry")
        self.tip_dia = _fan_geom_dict.get("tip_dia")
        self.htr = _fan_geom_dict.get("htr", 0.5)
        self.hub_dia = _fan_geom_dict.get("hub_dia",self.tip_dia*self.htr)
        self.sigma = _fan_geom_dict.get("sigma",1)
        
        _design_mean_line_dict = self._input_dict.get("turbomachinery",{}).get("design-mean-line")
        self.phi = _design_mean_line_dict.get("phi",0.7)
        self.eta_turbo = _design_mean_line_dict.get("eta_turbo",0.8)
        
        rcas = self.tip_dia/2 # radius of the casing
        rhub = self.hub_dia/2 # radius of the hub
        # _rm is the meanline radius in m at which mass flow above and below are equal
        self._rm = np.sqrt((rcas**2+rhub**2)/2)
        self._area1 = np.pi*(rcas**2-rhub**2) # area at the rotor
        
    
    def __str__(self):
        string = "Aero: {0}".format(self.name)
        string += "\n\tTip Diameter: {0} m".format(self.tip_dia)
        string += "\n\tHub Dia: {0} m".format(self.hub_dia)
        string += "\n\tFlow Coefficient: {0}".format(self.phi)
        return string
    
    # Getters and setters
    # -------------------
    
    def get_name(self):
        return self.name
    
    def get_area(self):
        return self._area1
    
    def get_tip_dia(self):
        return self.tip_dia
    
    def get_hub_dia(self):
        return self.hub_dia
    
    def set_tip_dia(self, tip_dia):
        self.tip_dia = tip_dia
    
    def set_hub_dia(self, hub_dia):
        self.hub_dia = hub_dia
    
    def set_phi(self, phi):
        self.phi = phi
    
    def set_sigma(self, sigma):
        self.sigma = sigma
    
    def set_etaFan(self, eta_turbo):
        self.eta_turbo = eta_turbo
    
    # Aerodynamic parameters for a given flight velocity and thrust
    # -------------------------------------------------------------
    
    def V1(self,rho, V0, thrust):
        """
        Calculates the flow velocity at the rotor for a given flight velocity and thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^2
        V0 : float
            Free stream velocity in m/s
        thrust : float
            Thrust in N

        Returns
        -------
        V1 : float
            Velocity at rotor in m/s

        """
        V1 = (V0+np.sqrt(V0**2+4*thrust/(rho*self._area1*self.sigma)))/(2/self.sigma)
        return V1
    
    def Vjet(self, rho, V0, thrust):
        """
        Calculate jet velocity for a given flight velocity and thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity in m/s
        thrust : float
            Thrust in N

        Returns
        -------
        Vjet : float
            Jet velocity

        """
        V1 = self.V1(rho, V0, thrust)
        Vjet = V1/self.sigma
        return Vjet
    
    def rotor_rpm(self, rho, V0, thrust):
        """
        Calculate rotor rpm required for a given flight velocity and thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity in m/s
        thrust : float
            Thrust in N

        Returns
        -------
        rpm : float
            Rotor RPM

        """
        V1 = self.V1(rho, V0, thrust)
        Umean = V1/self.phi
        rpm = Umean*60/(2*np.pi*self._rm)
        return rpm
    
    def shaftPower(self, rho, V0, thrust):
        """
        Calculates the power required for a given flight velocity and thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity in m/s
        thrust : float
            Thrust in N

        Returns
        -------
        float
            Shaft power in W

        """
        V1 = self.V1(rho, V0, thrust)
        Vjet = self.Vjet(rho, V0, thrust)
        power = 1/2*(Vjet**2*V1-V0**2*V1)*rho*self._area1
        return power/self.eta_turbo
    
    def torque(self, rho, V0, thrust):
        """
        Calculates the torque required for a given flight velocity and thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity in m/s
        thrust : float
            Thrust in N

        Returns
        -------
        float
            Shaft torque in N-m

        """
        rpm = self.rotor_rpm(rho, V0, thrust)
        shaftPower = self.shaftPower(rho, V0, thrust)
        return (shaftPower*60)/(2*np.pi*rpm)
    
    # Aerodynamic parameters for a given flight velocity and rpm
    # ----------------------------------------------------------
    
    def adv_ratio_J(self, V0, rpm):
        """
        Calculates advance ratio based on tip diameter

        Parameters
        ----------
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        J : float
            Advance Ratio based on tip diameter

        """
        n = rpm/60
        J = V0/(n*self.tip_dia)
        return J
    
    def adv_ratio_J1(self, V0, rpm):
        """
        Calculates advance ratio based on meanline radius

        Parameters
        ----------
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        J1 : float
            Advnace ratio based on meanline radius

        """
        J = self.adv_ratio_J(V0, rpm)
        J1 = J*self.tip_dia/(2*np.pi*self._rm)
        return J1
    
    def coef_thrust(self, V0, rpm):
        """
        Calculates coefficient of thrust

        Parameters
        ----------
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        CT : float
            Coefficient of Thrust

        """
        J1 = self.adv_ratio_J1(V0, rpm)
        CT1 = self.phi**2/self.sigma-self.phi*J1
        CT = CT1*self._area1/self.tip_dia**2*(2*np.pi*self._rm/self.tip_dia)**2
        return CT
    
    def coef_power(self, V0, rpm):
        """
        Calculates coefficient of power (shaft)

        Parameters
        ----------
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        CPshaft : float
            Ccoefficient of power (shaft)

        """
        J1 = self.adv_ratio_J1(V0, rpm)
        CP1 = 1/2*(self.phi**3/self.sigma**2-self.phi*J1**2)
        CP1shaft = CP1/self.eta_turbo
        CPshaft = CP1shaft*self._area1*(2*np.pi*self._rm)**3/self.tip_dia**5
        return CPshaft
    
    def coef_torque(self, V0, rpm):
        """
        Calculates coefficient of torque

        Parameters
        ----------
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        CTq : float
            Coefficient of torque

        """
        CPshaft = self.coef_power(V0, rpm)
        CTq = CPshaft/(2*np.pi)
        return CTq
    
    def thrust_from_rpm(self, rho, V0, rpm):
        """
        Calculate thrust from a given flight speed and rotor rpm

        Parameters
        ----------
        rho : float
            Density in kg/m^2
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        thrust : float
            Thrust in N

        """
        n = rpm/60
        CT = self.coef_thrust(V0, rpm)
        thrust = CT*rho*n**2*self.tip_dia**4
        return thrust
    
    def torque_from_rpm(self, rho, V0, rpm):
        """
        Calculates torque from a given flight speed and rotor rpm

        Parameters
        ----------
        rho : float
            Density in kg/m^2
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        torque : float
            Torque in N-m

        """
        n = rpm/60
        CTq = self.coef_torque(V0, rpm)
        torque = CTq*rho*n**2*self.tip_dia**5
        return torque
    
    def shaftPower_from_rpm(self, rho, V0, rpm):
        """
        Calculate shaft power from a given flight speed and rotor rpm

        Parameters
        ----------
        rho : float
            Density in kg/m^2
        V0 : float
            Free stream velocity in m/s
        rpm : float
            Rotor RPM

        Returns
        -------
        shaftPower : float
            Shaft power in W

        """
        n = rpm/60
        CPshaft = self.coef_power(V0, rpm)
        shaftPower = CPshaft*rho*n**3*self.tip_dia**5
        return shaftPower
    
    def Vjet_from_rpm(self, rpm):
        """
        Calculate jet velocity for a given rotor rpm

        Parameters
        ----------
        rpm : float
            Rotor rpm

        Returns
        -------
        TYPE
            Jet speed in m/s

        """
        return 2*np.pi*self._rm*self.phi*rpm/(60*self.sigma)
    
    # Other aerodynamic parameters
    # ----------------------------
    
    def thrust_from_Vjet(self, rho, V0, Vjet):
        """
        Calculate thrust for a given altitude, flight velocity and jet velocity

        Parameters
        ----------
        rho : float
            Density in kg/m^2
        V0 : float
            Free stream velocity in m/s
        Vjet : float
            Jet velocity in m/s

        Returns
        -------
        thrust : float
            Thrust in N

        """
        V1 = Vjet*self.sigma
        thrust = rho * self._area1 * V1 * (Vjet - V0)
        return thrust
    
    def psi_static(self):
        """
        Calculate stage loadin coefficinet, phsi, calculated from
        flow coefficient (phi) and sigma

        Returns
        -------
        float
            Stage loading coefficient

        """
        return self.phi**2/(2*self.sigma**2)
    
    def DH(self):
        """
        Calculates DH value

        Returns
        -------
        float
            DH value

        """
        psi = self.psi_static()
        return np.sqrt(self.phi**2+(1-psi)**2)/np.sqrt(self.phi**2+1)
    
    def eta_prop(self, V0):
        """
        Calculate propeller efficiency

        Parameters
        ----------
        V0 : float
            Free stream velocity

        Returns
        -------
        float
            Propeller efficiency

        """
        return 2*V0/(V0+self.Vjet_max)
    
    # Generate lists of two variables for plotting
    # --------------------------------------------
    
    def list_rpm_torque(self, rho, V0, rpm_max, num=40):
        """
        Generate lists to plot rpm vs torque

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        rpm_max : float
            Maximum rotor rpm
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        rpm_list : list
            List of rotor rpm
        torque_list : list
            List of torqe in N-m.

        """
        rpm_list = np.linspace(1000, rpm_max, num)
        torque_list = np.zeros(len(rpm_list))
        for i, rpm in enumerate(rpm_list):
            torque_list[i] = self.torque_from_rpm(rho, V0, rpm)
        return rpm_list, torque_list
    
    def list_rpm_thrust(self, rho, V0, rpm_max, num=40):
        """
        Generate lists to plot rpm vs thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        rpm_max : float
            Maximum rotor rpm
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        rpm_list : list
            List of rotor rpm
        thrust_list : list
            List of torqe in N-m.

        """
        rpm_list = np.linspace(1000, rpm_max, num)
        thrust_list = np.zeros(len(rpm_list))
        for i, rpm in enumerate(rpm_list):
            thrust_list[i] = self.thrust_from_rpm(rho, V0, rpm)
        return rpm_list, thrust_list
    
    def list_rpm_shaftPower(self, rho, V0, rpm_max, num=40):
        """
        Generate lists to plot rpm vs shaft power

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        rpm_max : float
            Maximum rotor rpm
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        rpm_list : list
            List of rotor rpm
        shaftPower_list : list
            List of shaft power in W.

        """
        rpm_list = np.linspace(1000, rpm_max, num)
        shaftPower_list = np.zeros(len(rpm_list))
        for i, rpm in enumerate(rpm_list):
            shaftPower_list[i] = self.shaftPower_from_rpm(rho, V0, rpm)
        return rpm_list, shaftPower_list
    
    def list_thrust_shaftPower(self, rho, V0, rpm_max, num=40):
        """
        Generate lists to plot thrust vs shaft power

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        rpm_max : float
            Maximum rotor rpm
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        thrust_list : list
            List of thrus in N
        shaftPower_list : list
            List of shaft power in W.

        """
        rpm_list = np.linspace(1000, rpm_max, num)
        thrust_list = np.zeros(len(rpm_list))
        shaftPower_list = np.zeros(len(rpm_list))
        for i, rpm in enumerate(rpm_list):
            thrust_list[i] = self.thrust_from_rpm(rho, V0, rpm)
            shaftPower_list[i] = self.shaftPower_from_rpm(rho, V0, rpm)
        return thrust_list, shaftPower_list
    
    # Plot two variables
    # ------------------
    
    def plot_rpm_thrust(self, rho, V0, num=40):
        """
        Plot rpm vs thrust

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        None.

        """
        rpm_list, thrust_list = self.list_rpm_thrust(rho, V0, num)
        fig, axes = plt.subplots()
        axes.plot(rpm_list, thrust_list)
        axes.set_xlabel("Rotor Speed [rpm]")
        axes.set_ylabel("Thrust [N]")
        axes.set_title("{}: RPM Vs Thrust at V0 = {:.1f} m/s".format(self.name, V0))
        fig.show()
    
    def plot_rpm_shaftPower(self, rho, V0, num=40):
        """
        Plot rpm vs shaft power

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        None.

        """
        rpm_list, shaftPower_list = self.list_rpm_shaftPower(rho, V0, num)
        fig, axes = plt.subplots()
        axes.plot(rpm_list, shaftPower_list)
        axes.set_xlabel("Rotor Speed [rpm]")
        axes.set_ylabel("Shaft Power [W]")
        axes.set_title("{}: RPM Vs Shaft Power at V0 = {:.1f} m/s".format(self.name, V0))
        fig.show()
    
    def plot_rpm_torque(self, rho, V0, num=40):
        """
        Plot rpm vs torque

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        None.

        """
        rpm_list, torque_list = self.list_rpm_torque(rho, V0, num)
        fig, axes = plt.subplots()
        axes.plot(rpm_list, torque_list)
        axes.set_xlabel("Rotor Speed [rpm]")
        axes.set_ylabel("Torque [N-m]")
        axes.set_title("{}: RPM Vs Torque at V0 = {:.1f} m/s".format(self.name, V0))
        fig.show()
    
    def plot_thrust_shaftPower(self, rho, V0, num=40):
        """
        Plot thrust vs shaft power

        Parameters
        ----------
        rho : float
            Density in kg/m^3
        V0 : float
            Free stream velocity
        num : int, optional
            Number of samples to generate. The default is 40. Must be non-negative

        Returns
        -------
        None.

        """
        thrust_list, shaftPower = self.list_thrust_shaftPower(rho, V0, num)
        fig, axes = plt.subplots()
        axes.plot(thrust_list, shaftPower)
        axes.set_xlabel("Thrust [N]")
        axes.set_ylabel("Shaft Power [W]")
        axes.set_title("{}: Thrust Vs Power at V0 = {:.1f} m/s".format(self.name, V0))
        fig.show()

# Testing
if __name__ == '__main__':
    # Define standard atmosphere
    from mission import standard_atmosphere
    std_atm = standard_atmosphere.StandardAtmosphere("SI")
    
    # Define Sycamore 160 aero
    aero_input = {
        "name": "Sycamore160",
        "turbomachinery": {
            "fan-geometry": {
                "tip_dia": 0.161,
    			"hub_dia": 0.057,
                "sigma": 1
            },
            "design-mean-line": {
                "phi": 0.581,
                "eta_turbo": 0.8
            }
        }
    }
    
    #aero_input = "design_param.json"
    sycamore160 = Aero1D(aero_input)
    print(sycamore160)
    
    # Some basic aero parameters
    altitude = 0
    rho = std_atm.rho(altitude)
    
    # Plot graphs at static condition
    V0 = 0
    rpm_max = 15000
    sycamore160.plot_rpm_shaftPower(rho, V0, rpm_max)
    sycamore160.plot_rpm_thrust(rho, V0, rpm_max)
    sycamore160.plot_rpm_torque(rho, V0, rpm_max)
    sycamore160.plot_thrust_shaftPower(rho, V0, rpm_max)
    
    V0 = 20
    sycamore160.plot_rpm_torque(rho, V0, rpm_max)