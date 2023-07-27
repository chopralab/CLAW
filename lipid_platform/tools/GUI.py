import os
import ipywidgets as widgets
import pandas as pd
import json



def folder_navigator():
    notebook_dir = os.getcwd()
    current_dir = os.path.join(notebook_dir, 'Projects')

    navigate_button = widgets.Button(description="Navigate")
    select_button = widgets.Button(description="Select this folder")
    select_current_button = widgets.Button(description="Select Current Folder")
    select = widgets.Select(options=[current_dir], rows=10)
    output = widgets.Output()

    Project_Folder = ""

    def on_navigate_button_clicked(b=None):
        nonlocal current_dir
        selected_option = select.value

        if selected_option == '..':
            current_dir = os.path.dirname(current_dir)
        else:
            current_dir = os.path.join(current_dir, selected_option)

        output.clear_output()
        dirs = next(os.walk(current_dir))[1]
        dirs.insert(0, '..') 
        select.options = dirs

        with output:
            print("Contents of {}: ".format(current_dir))
            print("\n".join(os.listdir(current_dir)))

    def on_select_button_clicked(b):
        nonlocal Project_Folder, current_dir
        selected_option = select.value

        if selected_option == '..':
            Project_Folder = os.path.dirname(current_dir)
        else:
            Project_Folder = os.path.join(current_dir, selected_option)
        output.clear_output()
        Project_Folder = Project_Folder + "/"

        with output:
            print("Selected Project Folder: ", Project_Folder)

        with open('Variable_Storage/folder_path.txt', 'w') as f:
            f.write(Project_Folder)

    def on_select_current_button_clicked(b):
        nonlocal Project_Folder, current_dir
        Project_Folder = current_dir

        output.clear_output()
        Project_Folder = Project_Folder + "/"
        with output:
            print("Selected Project Folder: ", Project_Folder)

        with open('Variable_Storage/folder_path.txt', 'w') as f:
            f.write(Project_Folder)

    navigate_button.on_click(on_navigate_button_clicked)
    select_button.on_click(on_select_button_clicked)
    select_current_button.on_click(on_select_current_button_clicked)

    display(navigate_button)
    display(select_button)
    display(select_current_button)
    display(select)
    display(output)

    # Programmatically display the contents of the initial directory
    on_navigate_button_clicked()



def load_project_folder():
    with open('Variable_Storage/folder_path.txt', 'r') as f:
        return f.read()


def assign_blank(unique_samples):
    dropdown = widgets.Dropdown(
        options=unique_samples,
        value=unique_samples[0],  # default value
        description='Samples',
        disabled=False,
    )
    display(dropdown)

    button = widgets.Button(description="Assign Blank")
    display(button)

    output = widgets.Output()
    display(output)

    blank_name = ""

    def on_button_clicked(b):
        nonlocal blank_name
        blank_name = dropdown.value

        with output:
            output.clear_output()
            print(f"Blank is: {blank_name}")

        # Save the blank_name to a file
        with open('Variable_Storage/blank_name.txt', 'w') as f:
            f.write(blank_name)

    button.on_click(on_button_clicked)

def load_blank_name():
    with open('Variable_Storage/blank_name.txt', 'r') as f:
        return f.read()




def filter_samples(labels_df):
    unique_samples = labels_df['Sample Name'].unique()

    multi_select = widgets.SelectMultiple(
        options=unique_samples,
        value=[unique_samples[0]],  # default value
        rows=len(unique_samples),
        description='Samples',
        disabled=False
    )
    display(multi_select)

    button = widgets.Button(description="Filter Samples")
    display(button)

    output = widgets.Output()
    display(output)

    labels_df2 = pd.DataFrame()

    def on_button_clicked(b):
        nonlocal labels_df2
        with output:
            output.clear_output()
        labels_df2 = labels_df[~labels_df['Sample Name'].isin(multi_select.value)]

        with output:
            display(labels_df2)

        # Save the DataFrame to a csv file
        labels_df2.to_csv('Variable_Storage/filtered_dataframe.csv')

    button.on_click(on_button_clicked)




def load_data_labels():
    # Specify the path to the csv file
    filtered_dataframe_file = 'Variable_Storage/filtered_dataframe.csv'
    
    # Check if the csv file exists and read it
    if os.path.isfile(filtered_dataframe_file):
        return pd.read_csv(filtered_dataframe_file, index_col=0)
    else:
        print("No saved DataFrame found.")
        return None
    


json_file_path = 'Variable_Storage/json_list_pairs.json'

def remove_empty_entries(json_list_pairs):
    cleaned_list_pairs = [
        [
            {key: value for key, value in pair_dict.items() if value} for pair_dict in pair
        ] for pair in json_list_pairs
    ]
    return cleaned_list_pairs

widgets_dict1 = {}
widgets_dict2 = {}

def display_pair_widgets(main_json):
    global widgets_dict1, widgets_dict2

    # Create a new global list each time widgets are displayed
    global json_list_pairs
    json_list_pairs = []

    # Check if file exists, if so delete it
    if os.path.exists(json_file_path):
        os.remove(json_file_path)

    widgets_dict1 = {key: widgets.SelectMultiple(options=value, description=key) for key, value in main_json.items()}
    widgets_dict2 = {key: widgets.SelectMultiple(options=value, description=key) for key, value in main_json.items()}
    
    for key in main_json.keys():
        display(widgets.HBox([widgets_dict1[key], widgets_dict2[key]]))

    def on_generate_clicked(b):
        new_json1 = {key: list(widget.value) for key, widget in widgets_dict1.items()}
        new_json2 = {key: list(widget.value) for key, widget in widgets_dict2.items()}
        pair = [new_json1, new_json2]
        json_list_pairs.append(pair)
        # Save data when Finish button is clicked
        save_data()

    def on_add_more_clicked(b):
        new_json1 = {key: list(widget.value) for key, widget in widgets_dict1.items()}
        new_json2 = {key: list(widget.value) for key, widget in widgets_dict2.items()}
        pair = [new_json1, new_json2]
        json_list_pairs.append(pair)
        
        for widget in widgets_dict1.values():
            widget.value = []
        for widget in widgets_dict2.values():
            widget.value = []

    generate_button = widgets.Button(description='Finish')
    generate_button.on_click(on_generate_clicked)

    add_more_button = widgets.Button(description='Add more JSON pairs')
    add_more_button.on_click(on_add_more_clicked)

    display(widgets.HBox([generate_button, add_more_button]))

def save_data():
    with open(json_file_path, 'w') as f:
        json.dump(json_list_pairs, f)

def load_data():
    with open(json_file_path, 'r') as f:
        return json.load(f)
    
def remove_empty_entries(json_list_pairs):
    cleaned_list_pairs = [
        [
            {key: value for key, value in pair_dict.items() if value} for pair_dict in pair
        ] for pair in json_list_pairs
    ]
    return cleaned_list_pairs


def get_unique_json_objects(json_list_pairs):
    json_set = set()
    for pair in json_list_pairs:
        for json_obj in pair:
            json_set.add(json.dumps(json_obj))
    
    json_list_singles = [json.loads(json_str) for json_str in json_set]
    return json_list_singles