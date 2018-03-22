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

Once you fill and press submit, it'll take some time and you'll be redirected to the results page where you can find our alignment along with [EBI's](https://www.ebi.ac.uk/) result for the same inputs and settings.

### Notes ###

- You must fill the gapopen and gapext parameters with negative values!
- If you don't fill theses parameters the default will be **-10 for gap** and **-0.5 for gapext**


##### django template from [@fasouto](https://github.com/fasouto/django-starter-template)