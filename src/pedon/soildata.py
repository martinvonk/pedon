
# import requests library
import requests
import pandas as pd

# send get request
r = requests.get(
    url = "https://www.soilphysics.wur.nl/soil.php",
    params = {"latitude": 52, "longitude": 5.5} )

# extract data in json format
data = r.json()

# further processing
pd.DataFrame(data["horizon"])
