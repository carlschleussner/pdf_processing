# PDF PROCESSING

"""
CLASS TO DERIVE REGIONAL AGGREGATION, PDF GENERATION AND KS TEST
FOLLOWING THE METHODOLOGY DEPLOYED
Schleussner, C.-F. et al.  ESD (2016)
http://www.earth-syst-dynam.net/7/327/2016/
by Carl Schleussner, Climate Analytics
carl.schleussner@climateanalytics.org
"""

# Requirements ESD processing class
    Input: 
        - Modelling and observed data for different resolution and time length
        - Needs to be NC dataframe with annual time axis 
    Output: 
        - Derive aggregation for reference period and target period(s)
        - SREX, Hemispheric and regionally resolved output 
        - Derive pdfs over periods and regional aggregation
        - Bootstrapping method for confidence intervals