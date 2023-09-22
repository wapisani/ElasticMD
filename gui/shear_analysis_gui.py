# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

Usage instructions:
    From the root of compchemkit, in a terminal or Anaconda prompt
    run the command "python -m pycct.gui.shear_analysis_gui"
    
Tested Python Versions:
    Python 3.8.13 (default, Mar 28 2022, 06:59:08) [MSC v.1916 64 bit (AMD64)]

"""

import os, io
from ..analysis import shear_functions as sf
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasAgg
import PySimpleGUI as sg

matplotlib.use("TkAgg")
sg.theme('Dark Blue 3')

def draw_figure(element, figure):
    """
    Draws the previously created "figure" in the supplied Image Element
    :param element: an Image Element
    :param figure: a Matplotlib figure
    :return: The figure canvas
    """

    plt.close('all')  # erases previously drawn plots
    canv = FigureCanvasAgg(figure)
    buf = io.BytesIO()
    canv.print_figure(buf, format='png')
    if buf is not None:
        buf.seek(0)
        element.update(data=buf.read())
        return canv
    else:
        return None


# Define the window layout
layout = [
    [sg.InputText(r"./ElasticMD/examples/PLA/Shear_Data/",font='Courier 18',size=(48,1),key="-directory-")],
    [sg.InputText(r"PLA_S2Sh1.log.lammps",font='Courier 18',size=(48,1),key="-logfile-")],
    [sg.Button("Load Data",font='Courier 18',key="-LOAD-",size=(12,1)),
    sg.Button("Fit Shear Modulus & Save Data",font='Courier 18',key="-fitShear-",size=(33,1),tooltip="Fitting may take several moments")],
    [sg.Text(text="Fitting parameters for G ",font='Courier 18',size=(45,1),key="-G_fit-")],
    [sg.Text(text="Number of breakpoints:",font='Courier 18',size=(41,1)),
     sg.InputText("2",font='Courier 18',size=(3,1),key="-GN-")],
    [sg.Text("Data range:",font='Courier 18',size=(31,1),key="-data_range_text-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-G_start-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-G_stop-")],
    [sg.Image(key="-G_CANVAS-")],
    [sg.Text(text="G=\nYield strain = \nYield Strength = ",font='Courier 18',size=(45,3),key="-G_result-")],
    
]

# Create the Window
window = sg.Window("Shear Analysis GUI", layout,finalize=True)

# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    print(event)
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break
    
    if event == '-LOAD-':
        # If you're loading a new file after completing an analysis,
        # clear the fitted figures
        if 'g_fig' in globals():
            plt.close(g_fig)
            window["-G_CANVAS-"].update()
        if 'fig_G' in globals():
            plt.close(fig_G)
            window["-G_CANVAS-"].update()
        
        log_file = values['-logfile-']
        os.chdir(os.path.normpath(values['-directory-']))
        strain_dir, primary_strain, primary_stress = sf.load_data(log_file)
        # Set parameter labels
        window['-G_fit-'].update(f"Fitting parameters for G{strain_dir}")
        window['-data_range_text-'].update(f"Data range: ({len(primary_strain)} points)")
        # Reset fitting parameters and results
        window['-GN-'].update("2")
        window['-G_start-'].update("0")
        window['-G_stop-'].update("0")
        window['-G_result-'].update(f"G{strain_dir}=\nYield strain = \nYield Strength = ")
        
        # Plot to canvases
        fig_G, ax_G = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_G.scatter(primary_strain,primary_stress,color="grey",s=20)
        ax_G.set_xlabel(f'True strain in {strain_dir}')
        ax_G.set_ylabel(f'True stress in {strain_dir}, MPa')
        draw_figure(window["-G_CANVAS-"], fig_G)
        
        
    if event == '-fitShear-':
        
        p_points = int(values["-GN-"])
        sh_start = int(values["-G_start-"])
        sh_end = int(values["-G_stop-"])
        shear_fit, shearMod, shearMod_stderror, yield_strain, yield_strength\
            = sf.compute_shear_modulus(primary_strain,primary_stress,p_points,data_range=(sh_start,sh_end))
            
        if shear_fit.get_results()['converged'] == True:
            # Erase the fitted figure when rerunning the analysis
            if 'g_fig' in globals():
                plt.close(g_fig)
                window["-G_CANVAS-"].update()
            else:
                plt.close(fig_G)
                window["-G_CANVAS-"].update()
                
            g_fig, g_ax = sf.plot_shear(shear_fit, strain_dir, shearMod, yield_strain, yield_strength)
            
            draw_figure(window["-G_CANVAS-"],g_fig)
            window["-G_result-"].update(f"G{strain_dir} = {shearMod:0.5f} +- {shearMod_stderror:0.5f} GPa\nYield strain = {yield_strain:0.4f}\nYield strength = {yield_strength:0.2f} MPa")
            
            window.refresh()
            
            sf.write_summaries(log_file,strain_dir,shearMod,yield_strain,\
                                yield_strength,shear_fit, shearMod_stderror)
            shear_filename = log_file.split('.log')[0]+f'_G{strain_dir}.pdf'
            g_fig.savefig(shear_filename,dpi=300)
            
            print(f"Data saved! Plot saved to {shear_filename}.")
    
window.close()
