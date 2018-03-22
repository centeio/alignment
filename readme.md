# SAAT

SUPER AWESOME ALIGNMENT TOOL

#### by [Carolina Centeio](https://github.com/centeio), [Eduarda Almeida](https://github.com/EAlmeida1) and [Higor Anjos](https://github.com/anjoshigor/)

## This is the web interface for the alignment python code done on the master brach of this repository ##

## Quickstart ##

Make sure to have python3 installed. Go to main folder and run:

    pip install requirements.txt

Once everything is installed you can run the development server: [http://localhost:8000/](http://localhost:8000/)

    python manage.py runserver

## How to use it ##

There's a pretty straightforward form you should fill with your sequences and the settings for the algorithms to work. **It's fasta format compatible!**

Once you fill and press "Alinhar", it'll take some time and you'll be redirected to the results page where you can find our alignment along with [EBI's](https://www.ebi.ac.uk/) result for the same inputs and settings.

### Usage Example

![Web Interface](https://im5.ezgif.com/tmp/ezgif-5-e48cb3bf8f.gif "Web Interface Usage")


##### django template from [@fasouto](https://github.com/fasouto/django-starter-template)