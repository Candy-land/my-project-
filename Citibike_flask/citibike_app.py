# -*- coding: utf-8 -*-
import gmplot
from flask import Flask, session, render_template, url_for,request
import matplotlib.pyplot
from flask_wtf import FlaskForm
from datetime import date
from wtforms.fields.html5 import DateField
from flaskext.mysql import MySQL
from wtforms import SelectField




app = Flask(__name__)

mysql = MySQL()
app.config['MYSQL_DATABASE_USER'] = 'root' 
app.config['MYSQL_DATABASE_PASSWORD'] = '190108' 
app.config['MYSQL_DATABASE_DB'] = 'citibike' 
app.config['MYSQL_DATABASE_HOST'] = 'localhost' 
mysql.init_app(app)
conn = mysql.connect() 
cursor = conn.cursor()

app.secret_key = '1234567890'

class DestinationsForm(FlaskForm):
    attributes = SelectField('Data Attributes', choices=[('Weekdays', 'Weekdays'), ('Weekends','Weekends'), ('Combined','Combined')])
"""
class AnalyticsForm(FlaskForm):
	attributes = SelectField('Data Attributes', choices=[('Agency', 'Agency'), ('Borough','Borough'), ('Complaint_Type', 'Complaint Type')])


class MapParamsForm(FlaskForm):
	dtfrom = DateField('DatePicker', format='%Y-%m-%d', default=date(2016,1,1)) 
	dtto = DateField('DatePicker', format='%Y-%m-%d', default=date(2016,1,2))
"""
def get_homepage_links():
	return [
			{"href": url_for('trips'), "label":"Trips Distribution in 2015"}, 
			{"href": url_for('gender'), "label":"Gender Distribution"},
			{"href": url_for('station'), "label":"Stations Growth"},
			{"href": url_for('destinations'), "label":"Most Popular Destinations"},
			{"href": url_for('age'), "label":"Age Distribution"},
			{"href": url_for('carto'), "label":"Close Look - Citibike Usage On 1/11/2016"},
			{"href": url_for('heatmap'), "label":"Heatmap"},
			]

@app.route("/")
def home():
    session["data_loaded"] = True
    return render_template('home.html',links=get_homepage_links())

@app.route("/station")#,methods=['GET','POST'])
def station():
	session["data_loaded"] = True
	return render_template('station.html')

@app.route("/trips")#,methods=['GET','POST'])
def trips():
	session["data_loaded"] = True
	return render_template('trips.html')

"""
@app.route("/destinations", methods=['GET','POST'])
def destinations():
	session["data_loaded"] = True
	return render_template('destinations.html',weekdayfile = 'weekday.html')
"""
@app.route("/destinations/", methods=['GET','POST'])#,methods=['GET','POST'])
def destinations():
    form = DestinationsForm()
    if form.validate_on_submit():
        #session["data_loaded"] = True
        if request.form['attributes'] == 'Weekdays':
            return render_template('weekdaybothpic.html')
        elif request.form['attributes'] == 'Weekends':
            return render_template('weekend_1bothpic.html')
        elif request.form['attributes'] == 'Combined':
            return render_template('weekday_weekend.html')
    return render_template('destinations.html', form=form)




@app.route("/growth")#,methods=['GET','POST'])
def growth():
	session["data_loaded"] = True
	return render_template('home.html')

@app.route("/carto")#,methods=['GET','POST'])
def carto():
	session["data_loaded"] = True
	return render_template('carto.html')

@app.route("/gender")#,methods=['GET','POST'])
def gender():
	session["data_loaded"] = True
	return render_template('gender.html')

@app.route("/heatmap")#,methods=['GET','POST'])
def heatmap():
	session["data_loaded"] = True
	return render_template('heatmap.html')

@app.route("/age")#,methods=['GET','POST'])
def age():
	session["data_loaded"] = True
	return render_template('age.html')

@app.route("/getstarted")#,methods=['GET','POST'])
def getstarted():
	session["data_loaded"] = True
	return render_template('getstarted.html')

@app.route("/references")#,methods=['GET','POST'])
def references():
	session["data_loaded"] = True
	return render_template('references.html')

#HTML PAGES
"""
@app.route("/weekday")#,methods=['GET','POST'])
def weekday():
	session["data_loaded"] = True
	return render_template('weekday.html')
"""
@app.route("/weekend")#,methods=['GET','POST'])
def weekend():
	session["data_loaded"] = True
	return render_template('weekend_1.html')

@app.route("/combined")#,methods=['GET','POST'])
def combined():
	session["data_loaded"] = True
	return render_template('Combine.html')

@app.route("/heat_nyc")#,methods=['GET','POST'])
def heat_nyc():
	session["data_loaded"] = True
	return render_template('heat_nyc.html')

@app.route("/weekday")#,methods=['GET','POST'])
def  weekday():
    session["data_loaded"] = True
    return render_template('weekday.html') 
    
@app.route("/weekend_1")#,methods=['GET','POST'])
def  weekend_1():
    session["data_loaded"] = True
    return render_template('weekend_1.html')


"""
def get_data(dtfrom,dtto):
	query = "select latitude, longitude from incidents where created_date >= '" + dtfrom + "' and created_date <= '" + dtto + "';" 
	cursor.execute(query) 
	return cursor.fetchall()

def get_df_data(): 
	import pandas
	query = "select unique_key, agency, complaint_type, borough from incidents;"
	cursor.execute(query)
	data = cursor.fetchall()
	df = pandas.DataFrame(data=list(data),columns=['Unique_key','Agency','Complaint Type','Borough']) 
	return df

@app.route("/map", methods=['GET','POST'])
def map():
	form = MapParamsForm()
	if form.validate_on_submit():
		dtfrom = form.dtfrom.data.strftime('%Y-%m-%d') 
		dtto = form.dtto.data.strftime('%Y-%m-%d') 
		coordinates = get_data(dtfrom, dtto)
		latitudes, longitudes = ([],[])
		if (len(coordinates)>0):
			for pair in coordinates: 
				latitudes.append(pair[0]) 
				longitudes.append(pair[1])
		gmap = gmplot.GoogleMapPlotter.from_geocode("New York",8) 
		gmap.heatmap(latitudes, longitudes) 
		gmap.draw('templates/mapoutput.html')
		return render_template('map.html', mapfile = 'mapoutput.html')
	return render_template('mapparams.html',form=form)

@app.route('/analytics/',methods=['GET','POST']) 
def analytics():
	form = AnalyticsForm()
	if form.validate_on_submit():
		import pandas
		df = get_df_data()
		column = request.form.get('attributes')
		group = df.groupby(column)
		ax = group.size().plot(kind='bar')
		fig = ax.get_figure() 
		fig.savefig('static/group_by_fig.png')
		return render_template('analyticsoutput.html')
	return render_template('analyticsparams.html', form=form)
"""

if __name__ == "__main__":
    app.run(debug=True)



