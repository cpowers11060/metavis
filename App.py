'''
Authors: Christopher Powers and Marvens Laporte
Institution: The University of Rhode Island


This code is the base for the app visualizing metabolic models and
allowing for the
'''


# import logging
# #===== START LOGGER =====
# logger = logging.getLogger(__name__)
# root_logger = logging.getLogger()
# root_logger.setLevel(logging.INFO)
# sh = logging.StreamHandler()
# formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# sh.setFormatter(formatter)
# root_logger.addHandler(sh)

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
from psamm.fluxanalysis import FluxBalanceProblem
from psamm.lpsolver import glpk, lp
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
        excluded_reactions=[]
        for rxn in nm.reactions:
            for cpd, v in rxn.equation.compounds:
                if (float(v)).is_integer()==False or float(v) > 10:
                    excluded_reactions.append(rxn.id)
        network = make_network_dict(nm, mm, subset=None, method='fpp',
                                                                element=el, excluded_reactions=excluded_reactions,
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
def build_network(nm, rxn_set, network, fba_dropdown):
        name={}
        formula={}
        for i in nm.compounds:
            name[i.id]=i.name
            formula[i.id]=i.formula
        nodes=[]
        edges=[]

        if not isinstance(fba_dropdown, list):
            objective=str(fba_dropdown)
            mm = nm.create_metabolic_model()

            # Set up the solver
            solver = glpk.Solver()
            # Define the flux balance
            problem = FluxBalanceProblem(mm, solver)
            problem.maximize(fba_dropdown)
            flux_carrying={}
            for i in nm.reactions:
                flux_carrying[i.id]=problem.get_flux(i.id)
        else:
            flux_carrying={}
            for i in nm.reactions:
                flux_carrying[i.id]=0


        for rxn in network[0]:
            if rxn.id in rxn_set:
                rxn_num=0
                for cpd in network[0][rxn][0]:
                    nodes.append({'data':{'id':str(cpd[0]),
                                                'label': name[str(cpd[0])[0:(str(cpd[0]).find('['))]],
                                                'formula': formula[str(cpd[0])[0:(str(cpd[0]).find('['))]]
                                                }})
                    nodes.append({'data':{'id':str(cpd[1]),
                                                'label': name[str(cpd[1])[0:(str(cpd[1]).find('['))]],
                                                'formula': formula[str(cpd[1])[0:(str(cpd[1]).find('['))]]
                                                }})
                    if  'pathways' in rxn.properties:
                        path = rxn.properties['pathways']
                    elif  'subsystem' in rxn.properties:
                        path = rxn.properties['subsystem']
                    else:
                        path = ['No pathway exists']
                    if rxn.id in edges:
                        edges.append({'data':{
                            'id':"".join([rxn.id, "_", str(rxn_num)]),
                            'source':str(cpd[0]),
                            'target':str(cpd[1]),
                            'label': "".join([rxn.name, "_", str(rxn_num)]),
                            'equation': str(rxn.properties["equation"]),
                            'pathways':path,
                            'flux':flux_carrying[rxn.id]
    #                            'equation':rxn.equation
                            }})
                        rxn_num+=1
                    else:
                        edges.append({'data':{
                            'id':"".join([rxn.id, "_", str(rxn_num)]),
                            'source':str(cpd[0]),
                            'target':str(cpd[1]),
                            'label': "".join([rxn.name, "_", str(rxn_num)]),
                            'equation': str(rxn.properties["equation"]),
                            'pathways':path,
                            'flux':flux_carrying[rxn.id]
    #                            'equation':rxn.equation
                            }})
                        rxn_num+=1

        return nodes, edges





# Generates all initial data for building the app
#nm, network = read_model("./models/E_rectale_MM/", "C")
#mr = ModelReader.reader_from_path("/Users/chrispowers/projects/ETH_Modelling/GEM-HS/model-hs.yaml")
#mr = ModelReader.reader_from_path("./models/E_rectale_MM/")
mr = ModelReader.reader_from_path("./models/iGEM_bin526_curated")
nm = mr.create_model()

pathway_list, rxn_set = get_pathway_list(nm, "All")
compounds_list = get_compounds_list(nm)
rxns=list(rxn_set)
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
                dbc.NavItem(
                        dbc.NavLink(
                                "Psamm documentation",
                                href="https://psamm.readthedocs.io",
                        )
                ),
                dbc.NavItem(
                        dbc.NavLink(
                                "Psamm publication",
                                href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004732",
                        )
                ),],

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
                        ##### Filter / Explore metabolic models
                        Use these filters to highlight reactions and compounds associated with different reactions
                        Try exploring different visualisation options.

                        The current visualization represents a draft model of Reinekea sp.
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
                                                            stylesheet=[
                                                                    {
                                                                            'selector': 'node',
                                                                            'style': {
                                                                                    'background-color': '#BFD7B5',
                                                                                    'label': 'data(label)'}},
                                                                    {
                                                                            "selector": "edge",
                                                                            "style": {
                                                                                    "width": 1,
                                                                                    "curve-style": "bezier"}},
                                                                    {
                                                                            "selector": "[flux != 0]",
                                                                            "style": {
                                                                                    "line-color": "blue"
                                                                            }
                                                                    }
                                                            ],
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
                                                        "Compound Search:", color="info", className="mr-1"),
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
                                                dbc.Badge(
                                                        "Flux Analysis:", color="info", className="mr-1"),
                                                dbc.FormGroup(
                                                        [
                                                                dcc.Dropdown(
                                                                        id="fba_dropdown",
                                                                        options=[
                                                                                {
                                                                                        "label": i,
                                                                                        "value": i,
                                                                                }
                                                                                for i in list(rxns)
                                                                        ],
                                                                        value=rxns,
                                                                        multi=False,
                                                                        style={"width": "100%"},
                                                                ),
                                                        ]
                                                ),
                                                dbc.Row(
                                                        [
                                                                dbc.Alert(
                                                                        id="fba-data",
                                                                        children="Select a reaction to see the flux here",
                                                                        color="secondary",
                                                                )
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
                        contents.append(
                                html.P(
                                        "Equation: "
                                        + str(data["equation"])
                                )
                        )
                        contents.append(
                                html.P(
                                        "Pathways: "
                                        + data["pathways"]
                                )
                        )
                        contents.append(
                                html.P(
                                        "Flux: "
                                        + str(data["flux"])
                                )
                        )


        return contents

@app.callback(
        Output("net", "elements"),
        [
                Input("pathways_dropdown", "value"),
                Input("element_dropdown", "value"),
                Input("compounds_dropdown", "value"),
                Input("fba_dropdown", "value"),
        ],
)
def filter_nodes(pathways_dropdown, element_dropdown, compounds_dropdown, fba_dropdown):

        #nm, network = read_model("/Users/chrispowers/projects/ETH_Modelling/GEM-HS/model-hs.yaml", element_dropdown)
        #nm, network = read_model("./models/E_rectale_MM/", element_dropdown)
        nm, network = read_model("./models/iGEM_bin526_curated", element_dropdown)
        pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)

        if isinstance(compounds_dropdown, str):
            rxn_list = []
            for rxn in network[0]:
                for cpd in network[0][rxn][0]:

                    if cpd[0].name == compounds_dropdown and rxn.id in rxn_set:
                        rxn_list.append(rxn.id)
        else:
            rxn_list = rxn_set
        nodes, edges = build_network(nm, rxn_list, network, fba_dropdown)
        elements=nodes+edges

        return elements


@app.callback(
        Output("fba-data", "children"), [Input("fba_dropdown", "value")]
)
def Run_FBA(fba_dropdown):
  if not isinstance(fba_dropdown, list):
      objective=str(fba_dropdown)
      # First read in the base model
      #mr = ModelReader.reader_from_path('../../ETH_Modelling/GEM-HS/E_rectale_MM/')
      #mr = ModelReader.reader_from_path("/Users/chrispowers/projects/ETH_Modelling/GEM-HS/model-hs.yaml")
      mr = ModelReader.reader_from_path("./models/iGEM_bin526_curated")
      nm = mr.create_model()
      mm = nm.create_metabolic_model()

      # Set up the solver
      solver = glpk.Solver()
      # Define the flux balance
      problem = FluxBalanceProblem(mm, solver)
      problem.maximize(fba_dropdown)



      return(str(problem.get_flux(fba_dropdown)))
  else:
      return("More than one reaction selected")


if __name__ == '__main__':
    app.run_server(debug=True)

#dashvis('Pyruvate metabolism' , 'C')



#anvil.server.wait_forever()
