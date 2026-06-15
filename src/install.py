import subprocess
import sys
import os
import importlib
import site

def install_and_import(install_name, import_name=None):
    if import_name is None:
        import_name = install_name
        
    try:
        # Try to import
        return importlib.import_module(import_name)
    except ImportError:
        print(f"{install_name} not found. Installing...")
        
        # Install the package
        subprocess.check_call([sys.executable, "-m", "pip", "install", install_name])
        
        # FORCE Python to look at the User Site Packages directory
        user_site = site.getusersitepackages()
        if user_site not in sys.path:
            sys.path.append(user_site)
        
        importlib.invalidate_caches()
        
        try:
            return importlib.import_module(import_name)
        except ImportError:
            # Specific fix for basemap which is a bit weird
            if install_name == 'basemap':
                return importlib.import_module('mpl_toolkits.basemap')
            raise

# --- Usage ---
# Standard packages
numpy = install_and_import('numpy')
pandas = install_and_import('pandas')
astropy = install_and_import('astropy')
sgp4 = install_and_import('sgp4')
xarray = install_and_import('xarray')
geopandas = install_and_import('geopandas')
shapely = install_and_import('shapely')
itur = install_and_import('itur')
numba = install_and_import('numba')

# SPECIAL CASE: Basemap
# It is installed as 'basemap' but lived inside 'mpl_toolkits'
try:
    from mpl_toolkits.basemap import Basemap
    print("Basemap imported successfully!")
except ImportError:
    install_and_import('basemap')
    from mpl_toolkits.basemap import Basemap
    print("Basemap installed and imported successfully!")

# Usage

install_and_import('basemap')
