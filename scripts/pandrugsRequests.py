import requests
from requests.auth import HTTPBasicAuth

sample_id = 'TCGA-05-4244-01A-01D-1105-08'


if __name__ == '__main__':

    # POST - computation
    url_post_computation = 'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/guest'
    # files = {'file': open(f'/home/cleon/extdata/input_pandrugs_{sample_id}.vcf', 'rb')}
    files = {'file': open(f'/home/cleon/extdata/tumor.vep_filtered.vcf', 'rb')}

    PARAMS = {'name': f'ComputationOf{sample_id}'}

    r = requests.post(url_post_computation, auth=HTTPBasicAuth('guest', 'guest'), params=PARAMS, files=files)
    print(r.text)

    # GET - variant analysis
    url_get_variant_analysis = f'https://www.pandrugs.org/pandrugs-backend/api/variantsanalysis/guest/{r.text}'

    r = requests.get(url=url_get_variant_analysis)
    data = r.json()
    print(data)

    '''

    # location given here
    location = "delhi technological university"

    # defining a params dict for the parameters to be sent to the API
    PARAMS = {'address': location}

    # sending get request and saving the response as response object
    r = requests.get(url=URL, params=PARAMS)

    # extracting data in json format
    data = r.json()

    # extracting latitude, longitude and formatted address
    # of the first matching location
    latitude = data['results'][0]['geometry']['location']['lat']
    longitude = data['results'][0]['geometry']['location']['lng']
    formatted_address = data['results'][0]['formatted_address']

    # printing the output
    print("Latitude:%s\nLongitude:%s\nFormatted Address:%s"
          % (latitude, longitude, formatted_address))'''
