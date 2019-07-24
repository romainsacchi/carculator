Validity tests
==============

Driving cycle, velocity and acceleration
----------------------------------------

There are eleven possible driving cycles to select from:

    * WLTC
    * WLTC 3.1
    * WLTC 3.2
    * WLTC 3.3
    * WLTC 3.4
    * CADC Urban
    * CADC Road
    * CADC Motorway
    * CADC Motorway 130
    * CADC
    * NEDC

Velocity, driving distance, driving time, and acceleration are calculated.

Manually, such parameters can be obtained the following way::
    
    import pandas as pd
    import numpy as np
    # Retrieve the driving cycle WLTC 3 from the UNECE
    driving_cycle = pd.read_excel('http://www.unece.org/fileadmin/DAM/trans/doc/2012/wp29grpe/WLTP-DHC-12-07e.xls',
                sheet_name='WLTC_class_3', skiprows=6, usecols=[2,4,5])

    # Calculate velocity (km/h -> m/s)
    velocity = driving_cycle['km/h'].values * 1000 / 3600

    # Retrieve driving distance (-> km)
    driving_distance = velocity.sum() * 1000

    # Retrieve driving time (-> s)
    driving_time = len(driving_cycle.values)

    # Retrieve acceleration by calculating the delta of velocity per time interval of 2 seconds
    acceleration = np.zeros_like(velocity)
    acceleration[1:-1] = (velocity[2:] - velocity[:-2])/2

Using `carculator`, these parameters can be obtained the following way::

    from coarse.energy_consumption import EnergyConsumptionModel
    ecm = EnergyConsumptionModel('WLTC')

    # Access the driving distance
    ecm.velocity.sum() * 1000

    # Access the driving time
    len(ecm.velocity)

    # Access the acceleration
    ecm.acceleration
    
Both approaches should return identical results::

    print(np.array_equal(velocity, ecm.velocity))
    print(driving_distance == ecm.velocity.sum()*1000)
    print(driving_time == len(ecm.velocity))
    print(np.array_equal(acceleration, ecm.acceleration))
    
    True
    True
    True
    True
    
And the acceleration returned by coarse should equal the values given by the UNECE::

    np.array_equal(np.around(ecm.acceleration,4),np.around(driving_cycle['m/s²'].values,4))
    
    True
    
Which can be also be verified visually::

    plt.plot(driving_cycle['m/s²'].values, label='UNECE')
    plt.plot(acceleration, label='Manually calculated')
    plt.plot(ecm.acceleration, label='coarse', alpha=0.6)
    plt.legend()
    plt.ylabel('m/s2')
    plt.xlabel('second')
    plt.savefig('comparison_driving_cycle.png')
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/comparison_driving_cycle.png
    :width: 400
    :alt: Alternative text
    
Car and components masses
-------------------------

This module returns the masses of the different vehicle components, given a vehicle size class, a year of manufacture and a powertrain techology. The relevant output is the `driving mass`, which will influence the energy required to overcome `rolling resistance`, the `drag`, but also the energy required to move the vehicle over a given distance -- `kinetic energy`.

Parameters such as total cargo mass, curb mass and driving mass, can be obtained the following way, for a 2017 battery electric SUV::

    cm.array.sel(size='SUV', powertrain='BEV', year=2017, parameter=['cargo mass','curb mass', 'driving mass']).values
    
    array([[  20.        ],
       [1719.56033224],
       [1874.56033224]])
       
One can check whether total cargo mass is indeed equal to cargo mass plus the product of the number of passengers and the average passenger weight::

    total_cargo, cargo, passengers, passengers_weight = cm.array.sel(size='SUV', powertrain='BEV', year=2017,
        parameter=['total cargo mass','cargo mass','average passengers', 'average passenger mass']).values
    print('Total cargo of {} kg, with a cargo mass of {} kg, and {} passengers of individual weight of {} kg.'.format(total_cargo[0], cargo[0], passengers[0], passengers_weight[0]))
    print(total_cargo == cargo+(passengers * passengers_weight))
    
    Total cargo of 155.0 kg, with a cargo mass of 20.0 kg, and 1.8 passengers of individual weight of 75.0 kg.
    [True]
    
However, most of the driving mass is explained by the curb mass:

    plt.pie(np.squeeze(cm.array.sel(size='SUV', powertrain='BEV', year=2017,
        parameter=['total cargo mass', 'curb mass']).values).tolist(), labels=['Total cargo mass', 'Curb mass'])
    plt.show()

.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/pie_total_mass.png
    :width: 400
    :alt: Alternative text
    
Here is a split between the components making up for the curb mass. One can see that, in the case of a battery electric SUV, most of the weight comes from the glider as well as the battery cells. On an equivalent diesel powertrain, the mass of the gliber base is comparatively more important::

    l_param=["fuel mass","charger mass","converter mass","glider base mass","inverter mass","power distribution unit mass",
            "combustion engine mass","electric engine mass","powertrain mass","fuel cell stack mass",
            "fuel cell ancillary BoP mass","fuel cell essential BoP mass","battery cell mass","battery BoP mass","fuel tank mass"]


    colors = ['yellowgreen','red','gold','lightskyblue','white','lightcoral','blue','pink', 'darkgreen','yellow','grey','violet','magenta','cyan', 'green']

    BEV_mass = np.squeeze(cm.array.sel(size='SUV', powertrain='BEV', year=2017,
            parameter=l_param).values)

    percent = 100.*BEV_mass/BEV_mass.sum()

    f = plt.figure(figsize=(15,10))

    ax = f.add_subplot(121)

    patches, texts = ax.pie(BEV_mass, colors=colors, startangle=90, radius=1.2)
    ax.set_title('BEV SUV')
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(l_param, percent)]

    sort_legend = True
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, BEV_mass),
                                              key=lambda x: x[2],
                                              reverse=True))

    ax.legend(patches, labels, loc='upper left', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)


    ICEV_d_mass = np.squeeze(cm.array.sel(size='SUV', powertrain='ICEV-d', year=2017,
            parameter=l_param).values)
    percent = 100.*ICEV_d_mass/ICEV_d_mass.sum()

    ax2 = f.add_subplot(122)

    patches, texts = ax2.pie(ICEV_d_mass, colors=colors, startangle=90, radius=1.2)
    ax2.set_title('ICE-d SUV')
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(l_param, percent)]

    sort_legend = True
    if sort_legend:
        patches, labels, dummy =  zip(*sorted(zip(patches, labels, ICEV_d_mass),
                                              key=lambda x: x[2],
                                              reverse=True))

    ax2.legend(patches, labels, loc='upper left', bbox_to_anchor=(-0.1, 1.),
               fontsize=8)

    plt.subplots_adjust(wspace=1)
    plt.show()
  
.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/pie_mass_components.png
    :width: 400
    :alt: Alternative text
    

Model vs. car manufacturers' data
*********************************

The `curb mass` returned by Carculator is plotted against manufacturers' data, per vehicle size class and powertrain technology. To do so, we use the car database Car2db (https://car2db.com/) and load all car trims produced after 2013 (21,383 vehicles)::

    import mysql.connector
    import matplotlib.pyplot as plt


    cnx = mysql.connector.connect(user='root', password='****',
                                  host='localhost',
                                  database='cardb')
    cursor = cnx.cursor()

    year = 2013

    query_ID = ("SELECT DISTINCT ct.id_car_trim "
             "FROM car_make c "
                 "INNER JOIN car_model cm ON cm.id_car_make = c.id_car_make "
                 "INNER JOIN car_generation cg ON cg.id_car_model = cm.id_car_model "
                 "INNER JOIN car_trim ct ON ct.id_car_model = cm.id_car_model "
                 "INNER JOIN car_specification_value csv ON csv.id_car_trim = ct.id_car_trim "
                 "INNER JOIN car_specification  cs ON cs.id_car_specification = csv.id_car_specification "
                 "WHERE ct.start_production_year > %s OR cg.year_begin > %s"    
            )

    cursor.execute(query_ID, (year, year))
    l_ID = [l[0] for l in cursor.fetchall()]

    cnx.close()

    cnx = mysql.connector.connect(user='root', password='*****',
                                  host='localhost',
                                  database='cardb')

    cursor = cnx.cursor()

    query = ("SELECT ct.id_car_trim, c.name, cm.name, ct.name, cs.name, csv.value, ct.start_production_year "
             "FROM car_make c "
                 "INNER JOIN car_model cm ON cm.id_car_make = c.id_car_make "
                 "INNER JOIN car_trim ct ON ct.id_car_model = cm.id_car_model "
                 "INNER JOIN car_specification_value csv ON csv.id_car_trim = ct.id_car_trim "
                 "INNER JOIN car_specification  cs ON cs.id_car_specification = csv.id_car_specification "
            "WHERE  ct.id_car_trim IN (" + ','.join(map(str, l_ID)) + ")"
            )
    cursor.execute(query)

    df = pd.DataFrame( [[ij for ij in i] for i in cursor.fetchall()])
    df.rename(columns={0:'ID', 1: 'Make name', 2: 'Model name', 3:'Trim', 4:'Spec', 5:'Spec value', 6:'Year'}, inplace=True)

    cnx.close()

    df = df.pivot_table(index=['ID', 'Make name', 'Model name', 'Trim'], columns='Spec', values='Spec value', aggfunc='first')
    df['Volume'] = (df['Height'] * df['Width'] * df['Length']) / 1e9
    df['Footprint'] = ((df['Front track'] + df['Rear track']) * df['Wheelbase'])/2/1e6
    df = df[(~df['Footprint'].isnull())&(~df['Curb weight'].isnull())]
    
    def find_class(row):
        if row['Body type'] == 'Pickup' or row['Body type'] == 'Minivan': return 'Van'
        if row['Body type'] == 'Crossover': return 'SUV'
        if row['Footprint'] <= 3.4 and row['Curb weight'] <= 1050: return 'Mini'
        if row['Footprint'] > 3.4 and row['Footprint'] <= 3.8 and row['Curb weight'] > 900 and row['Curb weight'] <= 1250: return 'Small'
        if row['Footprint'] > 3.8 and row['Footprint'] <= 4.3 and row['Curb weight'] > 1250 and row['Curb weight'] <= 1500: return 'Lower medium'
        if row['Footprint'] > 4.1 and row['Footprint'] <= 4.4 and row['Curb weight'] > 1450 and row['Curb weight'] <= 1750: return 'Medium'
        if row['Footprint'] > 4.4 and row['Curb weight'] > 1450 : return 'Large'

    df.loc[:,'coarse class']=''
    df.loc[:,['coarse class']] = df.apply(find_class, axis=1)
    
    plt.style.use('seaborn')
    d_pt = {'Diesel':'ICEV-d',
           'Electric':'BEV',
           'Gas':'ICEV-g',
           'Gasoline':'ICEV-p',
           'Gasoline, Gas':'ICEV-g',
           'Hybrid':'HEV-p',
           'Diesel, Hybrid': 'HEV-p',
           'Gasoline, Electric': 'HEV-p',
           'Fuel cell': 'FCEV'}

    df_graph = df
    df_graph = df_graph.append({'Make name':'Hyundai','Model name':'ix35','Engine type': 'Fuel cell', 'coarse class':'SUV', 'Curb weight': 1846}, ignore_index=True)
    df_graph = df_graph.append({'Make name':'Toyota','Model name':'Mirai','Engine type': 'Fuel cell', 'coarse class':'Medium', 'Curb weight': 1850}, ignore_index=True)
    df_graph = df_graph.append({'Make name':'Honda','Model name':'Clarity','Engine type': 'Fuel cell', 'coarse class':'Large', 'Curb weight': 1876}, ignore_index=True)
    df_graph = df_graph.replace({"Engine type": d_pt})
    grouped = df_graph.groupby(['Engine type'])
    rowlength = grouped.ngroups                     # fix up if odd number of groups
    plt.close('all')
    fig, axs = plt.subplots(figsize=(16,8), 
                            nrows=2, ncols=3,     # fix as above
                            gridspec_kw=dict(hspace=0.4),
                           sharey=True) # Much control of gridspec

    targets = zip(grouped.groups.keys(), axs.flatten())


    for i, (key, ax) in enumerate(targets):
        order =  ['Mini', 'Small', 'Lower medium', 'Medium', 'Large','SUV', 'Van']

        sns.boxplot(y='Curb weight', x='coarse class', 
                     data=grouped.get_group(key), 
                     width=0.5, showfliers=False,
                     color='green', ax=ax, zorder=5, order =order)

        np.squeeze(cm.array.sel(size=order, powertrain=key, year=2017,
            parameter='curb mass')).to_dataframe(name='weight').plot(y='weight',ax=ax, kind='line', linestyle='None', marker='o',
                                                                    markerfacecolor='red', markersize=10, zorder=10)

        ax2 = grouped.get_group(key).groupby('coarse class').size().reindex(index = [o for o in order]).plot.bar(secondary_y = True,
                                                                            ax = ax, alpha = 0.2, zorder=0)
        ax.set_title(key + ' (' + str(grouped.get_group(key).groupby('coarse class').size().sum()) + ' datasets)')

        if (i ==2 or i==5):
            ax2.set_ylabel('Number of datasets', rotation=270, labelpad=20)


        ax.set_ylim(0,2650)
        ax2.grid(False)
        ax2.set_ylim(0,np.max(grouped.get_group(key).groupby('coarse class').size()*1.2))

        for tick in ax.get_xticklabels():
            tick.set_rotation(45)

    red_patch = mpatches.Patch(color='red', label='Carculator value')
    blue_patch = mpatches.Patch(color='steelblue', alpha=.2, label='number of vehicles in Car2db')
    bar_patch =  mpatches.Patch(color='green', label='Car2db values')
    plt.legend(handles=[red_patch, bar_patch, blue_patch],bbox_to_anchor=(-1.8,-.4))

    for a in fig.axes:
        a.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=True,
        top=False,
        labelbottom=True,
        rotation=45)    # labels along the bottom edge are on
    plt.tight_layout()
    plt.show()  
    
.. image:: https://github.com/romainsacchi/coarse/raw/master/docs/mass_comparison.png
    :width: 400
    :alt: Alternative text
    
    

