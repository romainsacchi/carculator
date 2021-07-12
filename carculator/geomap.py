from wurst.geo import geomatcher

from . import DATA_DIR

REGION_MAPPING_FILEPATH = DATA_DIR / "regionmappingH12.csv"


def get_IAM_geomatcher():
    """
    Geographical boundaries for IMAGE regions are initally included in geomatcher.
    However, they are not properly labelled.

    """

    d_image_regions = {
        "BRA": "Brazil",
        "CAN": "Canada",
        "CEU": "Central Europe",
        "CHN": "China Region",
        "EAF": "Eastern Africa",
        "INDIA": "India",
        "INDO": "Indonesia Region",
        "JAP": "Japan",
        "KOR": "Korea Region",
        "ME": "Middle east",
        "MEX": "Mexico",
        "NAF": "Northern Africa",
        "OCE": "Oceania",
        "RCAM": "Central America",
        "RSAF": "Rest of Southern Africa",
        "RSAM": "Rest of South America",
        "RSAS": "Rest of South Asia",
        "RUS": "Russia Region",
        "SAF": "South Africa",
        "SEAS": "South Asia",
        "STAN": "Central Asia",
        "TUR": "Turkey",
        "UKR": "Ukraine region",
        "USA": "USA",
        "WAF": "Western Africa",
        "WEU": "Western Europe",
    }

    d_map = {("IMAGE", v): ("IMAGE", k) for k, v in d_image_regions.items()}

    new_def = dict()

    for k, v in geomatcher.items():
        if isinstance(k, tuple):
            if k[0] == "IMAGE" and k[1] in list(d_image_regions.values()):
                new_def[d_map[k]] = v

    geo = geomatcher

    for k in list(geomatcher.keys()):
        if k[0] == "IMAGE" and k[1] in list(d_image_regions.values()):
            geomatcher.pop(k)

    geo.update(new_def)

    with open(REGION_MAPPING_FILEPATH) as f:
        f.readline()
        csv_list = [[val.strip() for val in r.split(";")] for r in f.readlines()]
        l = [(x[1], x[2]) for x in csv_list]

    # List of countries not found
    countries_not_found = ["CC", "CX", "GG", "JE", "BL"]

    rmnd_to_iso = {}
    iso_to_rmnd = {}

    # Build a dictionary that maps region names (used by REMIND) to ISO country codes
    # And a reverse dictionary that maps ISO country codes to region names
    for ISO, region in l:
        if ISO not in countries_not_found:
            try:
                rmnd_to_iso[region].append(ISO)
            except KeyError:
                rmnd_to_iso[region] = [ISO]

            iso_to_rmnd[region] = ISO

    geo.add_definitions(rmnd_to_iso, "REMIND")

    return geo


class Geomap:
    """
    Map ecoinvent locations to IAM regions and vice-versa.
    """

    def __init__(self):

        self.geo = get_IAM_geomatcher()

    def iam_to_ecoinvent_location(self, location, contained=False):
        """
        Find the corresponding ecoinvent region given an IAM region.

        :param location: name of a IAM region
        :type location: str
        :return: name of an ecoinvent region
        :rtype: str
        """

        if location == "World":
            return ["GLO"]

        ecoinvent_locations = []

        searchfunc = self.geo.contained if contained else self.geo.intersects

        for iam in ("REMIND", "IMAGE"):
            loc = (iam, location)

            try:
                searchfunc(loc)
                for r in searchfunc(loc):
                    if not isinstance(r, tuple):
                        ecoinvent_locations.append(r)

            except KeyError:
                pass

        if len(ecoinvent_locations) == 0:
            print("Can't find location {} using the geomatcher.".format(location))

        return ecoinvent_locations
