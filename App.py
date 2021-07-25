'''
Authors: Christopher Powers and Marvens Laporte
Institution: The University of Rhode Island


This code is the base for the app visualizing metabolic models and
allowing for the
'''


import logging
#===== START LOGGER =====
logger = logging.getLogger(__name__)
root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
sh = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
sh.setFormatter(formatter)
root_logger.addHandler(sh)

from collections import defaultdict
from psamm.datasource.native import NativeModel, ModelReader, ModelWriter
from psamm.expression import boolean
from psamm.lpsolver import generic
from psamm.fluxanalysis import FluxBalanceProblem
from psamm.lpsolver import glpk, lp
from psamm.graph import make_network_dict, write_network_dict
from psamm.datasource.reaction import Reaction, Compound
import plotly.express as px
import sys

# From Bokeh/NetworkX


# From Visdxx
import json
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import dash_cytoscape as cyto


el='C'
pathway='All'

# Generates the model and the findprimarypairs network
def read_model(model_path, el):
        if isinstance(el, list) or el == None:
            el = "C"
        # First read in the base model
        mr = ModelReader.reader_from_path(model_path)
        nm = mr.create_model()
        mm = nm.create_metabolic_model()

        network = make_network_dict(nm, mm, subset=None, method='fpp',
                                                                element=el, excluded_reactions=[],
                                                                reaction_dict={}, analysis=None)
        return nm, network


# Builds a reaction set from a model based on a pathway of interest
def get_pathway_list(nm, pathway):
        pathway_list=set()
        if isinstance(pathway, list) or pathway == None:
            pathway = "All"
        if pathway != "All":
            rxn_set=set()
            for i in nm.reactions:
                if 'pathways' in i.properties:
                    for j in i.properties['pathways']:
                        pathway_list.add(j)
                    if pathway in i.properties['pathways']:
                        rxn_set.add(i.id)
                elif  'subsystem' in i.properties:
                    pathway_list.add(i.properties['subsystem'])
                    if pathway in i.properties['subsystem']:
                        rxn_set.add(i.id)

        else:
            rxn_set=set()
            for i in nm.reactions:
                if  'pathways' in i.properties:
                    for j in i.properties['pathways']:
                        pathway_list.add(j)
                elif  'subsystem' in i.properties:
                    pathway_list.add(i.properties['subsystem'])

                rxn_set.add(i.id)
        pathway_list=["All"] + list(pathway_list)
        return pathway_list, rxn_set

def get_compounds_list(nm):
    compounds_list=[]
    for i in nm.compounds:
        compounds_list.append(i.id)
    return compounds_list
# Useful function to build the subset network. Returns nodes and edges
# from the network associated with the rxn_set of interest
def build_network(nm, rxn_set, network):
        name={}
        formula={}
        for i in nm.compounds:
            name[i.id]=i.name
            formula[i.id]=i.formula
        nodes=[]
        edges=[]
        for rxn in network[0]:
            if rxn.id in rxn_set:
                for cpd in network[0][rxn][0]:
                    nodes.append({'data':{'id':str(cpd[0]),
                                                'label': name[str(cpd[0])[0:-3]],
                                                'formula': formula[str(cpd[0])[0:-3]]
                                                }})
                    nodes.append({'data':{'id':str(cpd[1]),
                                                'label': name[str(cpd[1])[0:-3]],
                                                'formula': formula[str(cpd[1])[0:-3]]
                                                }})
                    if  'pathways' in rxn.properties:
                        path = rxn.properties['pathways']
                    elif  'subsystem' in rxn.properties:
                        path = rxn.properties['subsystem']
                    else:
                        path = ['No pathway exists']
                    edges.append({'data':{
                            'id':rxn.id,
                            'source':str(cpd[0]),
                            'target':str(cpd[1]),
                            'label': rxn.name,
                            'pathways':path
#                            'equation':rxn.equation
                            }})
        return nodes, edges



# Generates all initial data for building the app
nm, network = read_model("./models/iGEM_bin526_curated/", "C")
pathway_list, rxn_set = get_pathway_list(nm, "All")
compounds_list = get_compounds_list(nm)
# nodes, edges = build_network(nm, rxn_set, network)
# initialize an empty list. the full is messy
nodes, edges = [], []


# Initialize the app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

# Set styles for the app
styles = {
        'pre': {
                'border': 'thin lightgrey solid',
                'overflowX': 'scroll'
        }
}
col_swatch = px.colors.qualitative.Dark24
default_stylesheet = [
        {
                'selector': 'node',
                'style': {
                        'background-color': '#BFD7B5',
                        'label': 'data(label)'}},
        {
                "selector": "edge",
                "style": {
                        "width": 1,
                        "curve-style": "bezier"}}
]

# Set the navigation bar for the app
navbar = dbc.NavbarSimple(
        children=[
                dbc.NavItem(
                        dbc.NavLink(
                                "Source Code",
                                href="https://github.com/cpowers11060/metavis",
                        )
                ),
        ],
        brand="Psamm web-based visualization of metabolic models",
        brand_href="#",
        color="dark",
        dark=True,
)

# Define the body of the app
body_layout = dbc.Container(
        [
                dbc.Row(
                        [
                                dcc.Markdown(
                                        """
                        -----
                        ##### Filter / Explore metabolic model of _Reinekea_ MAG from NB
                        Use these filters to highlight reactions and compounds associated with different reactions
                        Try exploring different visualisation options.
                        -----
                        """
                                ),
                        ]
                ),
                dbc.Row(
                        [
                                dbc.Col(
                                        [
                                                dbc.Row(
                                                        [
                                                            cyto.Cytoscape(id = 'net',
                                                            layout={'name':'cose'},
                                                            style={'width': '500px', 'height': '500px'},
                                                            elements=nodes+edges,
                                                            stylesheet=default_stylesheet,
                                                            minZoom=0.06
                                                            )
                                                        ]
                                                ),

                                        ],
                                        sm=12,
                                        md=8,
                                ),
                                dbc.Col(
                                        [

                                                dbc.Row(
                                                        [
                                                                dbc.Alert(
                                                                        id="node-data",
                                                                        children="Click on a node to see its details here",
                                                                        color="secondary",
                                                                )
                                                        ]
                                                ),
                                                dbc.Row(
                                                        [
                                                                dbc.Alert(
                                                                        id="edge-data",
                                                                        children="Click on an edge to see its details here",
                                                                        color="secondary",
                                                                )
                                                        ]
                                                ),
                                                dbc.Badge(
                                                        "Pathways:", color="info", className="mr-1"
                                                ),
                                                dbc.FormGroup(
                                                        [
                                                                dcc.Dropdown(
                                                                        id="pathways_dropdown",
                                                                        options=[
                                                                                {
                                                                                        "label": i,
                                                                                        "value": i,
                                                                                }
                                                                                for i in list(pathway_list)
                                                                        ],
                                                                        value=pathway_list,
                                                                        multi=False,
                                                                        style={"width": "100%"},
                                                                ),
                                                        ]
                                                ),
                                                dbc.Badge(
                                                        "Element Transfer Networks:", color="info", className="mr-1"
                                                ),
                                                dbc.FormGroup(
                                                        [
                                                                dcc.Dropdown(
                                                                        id="element_dropdown",
                                                                        options=[
                                                                                {
                                                                                        "label": i,
                                                                                        "value": i,
                                                                                }
                                                                                for i in ["C","N","S","P"]
                                                                        ],
                                                                        value=["C","N","S","P"],
                                                                        multi=False,
                                                                        style={"width": "100%"},
                                                                ),
                                                        ]
                                                ),
                                                dbc.Badge(
                                                        "Compound Networks:", color="info", className="mr-1"),
                                                dbc.FormGroup(
                                                        [
                                                                dcc.Dropdown(
                                                                        id="compounds_dropdown",
                                                                        options=[
                                                                                {
                                                                                        "label": i,
                                                                                        "value": i,
                                                                                }
                                                                                for i in list(compounds_list)
                                                                        ],
                                                                        value=compounds_list,
                                                                        multi=False,
                                                                        style={"width": "100%"},
                                                                ),
                                                        ]
                                                ),


                                                dbc.Badge(color="info", className="mr-1"),
                                                dbc.FormGroup(
                                                        [
                                                                dbc.Container(
                                                                        [
                                                                                dbc.Checkbox(
                                                                                        id="show_edges_radio",
                                                                                        className="form-check-input",
                                                                                        checked=True,
                                                                                ),
                                                                                dbc.Label(
                                                                                        "Show Reactions",
                                                                                        html_for="show_edges_radio",
                                                                                        className="form-check-label",
                                                                                        style={
                                                                                                "color": "DarkSlateGray",
                                                                                                "fontSize": 12,
                                                                                        },
                                                                                ),
                                                                        ]
                                                                )
                                                        ]
                                                ),
                                        ],
                                        sm=12,
                                        md=4,
                                ),
                        ]
                ),
                dbc.Row(
                        [
                                dcc.Markdown(
                                        """
                        \* Data analysis carried out for demonstration of data visualisation purposes only.
                        """
                                )
                        ],
                        style={"fontSize": 11, "color": "gray"},
                ),
        ],
        style={"marginTop": 20},
)

app.layout = html.Div([navbar, body_layout])

@app.callback(
        Output("node-data", "children"), [Input("net", "selectedNodeData")]
)
def display_nodedata(datalist):
        contents = "Click on a node to see its details here"
        if datalist is not None:
                if len(datalist) > 0:
                        data = datalist[-1]
                        contents = []
                        contents.append(
                                html.H5(
                                                "ID: " + data["id"].title()
                                )
                        )
                        contents.append(
                                html.P(
                                        "Name: " + data["label"]
                                )
                        )
                        contents.append(
                                html.P(
                                        "Formula: "
                                        + str(data["formula"])
                                )
                        )

        return contents

@app.callback(
        Output("edge-data", "children"), [Input("net", "selectedEdgeData")]
)
def display_nodedata(datalist):
        contents = "Click on an edge to see its details here"
        if datalist is not None:
                if len(datalist) > 0:
                        data = datalist[-1]
                        contents = []
                        contents.append(html.H5("ID: " + data["id"].title()))
                        contents.append(
                                html.P(
                                        "Name: "
                                        + data["label"]
                                )
                        )
                        # contents.append(
                        #         html.P(
                        #                 "Equation: "
                        #                 + str(data["equation"])
                        #         )
                        # )
                        contents.append(
                                html.P(
                                        "Pathways: "
                                        + data["pathways"]
                                )
                        )

        return contents

@app.callback(
        Output("net", "elements"),
        [
                Input("pathways_dropdown", "value"),
                Input("element_dropdown", "value"),Input("compounds_dropdown", "value"),
        ],
)
def filter_nodes(pathways_dropdown, element_dropdown, compounds_dropdown):

        nm, network = read_model("./models/iGEM_bin526_curated/", element_dropdown)
        pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)
        print(type(compounds_dropdown))
        if isinstance(compounds_dropdown, str):
            rxn_list = []
            for rxn in network[0]:
                for cpd in network[0][rxn][0]:
                    print(cpd[0].name, compounds_dropdown)
                    if cpd[0].name == compounds_dropdown and rxn.id in rxn_set:
                        rxn_list.append(rxn.id)
        else:
            rxn_list = rxn_set
        nodes, edges = build_network(nm, rxn_list, network)
        elements=nodes+edges

        return elements

# @app.callback(
#         Output("compunds net", "compounds"),
#         [
#                 Input("compounds_dropdown", "value"),
#                 Input("element_dropdown", "value"),
#         ],
# )
# def filter_compounds(compounds_dropdown, element_dropdown):
#
#         nm, network = read_model("./models/iGEM_bin526_curated/", element_dropdown)
#         rxn_list = []
#         for rxn in network[0]:
#             for cpd in network[0][rxn][0]:
#                 if cpd[0].name == compounds_dropdown:
#                     rxn_list.append(rxn.id)
#         nodes, edges = build_network(nm, rxn_list, network)
#         compounds=nodes+edges
#
#         return compounds

if __name__ == '__main__':
    app.run_server(debug=True)

#dashvis('Pyruvate metabolism' , 'C')



#anvil.server.wait_forever()
