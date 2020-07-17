from lib.aerodynamics import Aerodynamics

aerodynamics = Aerodynamics()
aerodynamics.setAerodynmamicsData()
withem = aerodynamics.getWithem()
for key in withem:
    print('key = ', key, '________witham = ', withem[key])