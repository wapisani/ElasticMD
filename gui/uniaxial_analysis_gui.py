# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

Usage instructions:
    From the root of compchemkit, in a terminal or Anaconda prompt
    run the command "python -m ElasticMD.gui.uniaxial_analysis_gui"
    
Tested Python Versions:
    Python 3.8.13 (default, Mar 28 2022, 06:59:08) [MSC v.1916 64 bit (AMD64)]
"""

import os, io
from ..analysis import uniaxial_functions as uf
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
    [sg.InputText(r"C:\Users\RDEL1WAP\Documents\Research\GSL\Aerogels\PU",font='Courier 18',size=(45,1),key="-directory-"),
     sg.InputText(r"n3300_300K_2x2.Comp2.log.lammps",font='Courier 18',size=(45,1),key="-logfile-"),
     sg.Button("Load Data",font='Courier 18',key="-LOAD-"),
     sg.Button("Fit All & Save Data",font='Courier 18',key="-FITALL-"),
     sg.Button("Save Data",font='Courier 18',key="-SAVE-")],
    [sg.Text(text="Fitting parameters for E ",font='Courier 18',size=(45,1),key="-E_fit-"),
     sg.Text(text="Fitting parameters for Nu  ",font='Courier 18',size=(45,1),key="-Nu1_fit-"),
     sg.Text(text="Fitting parameters for Nu  ",font='Courier 18',size=(43,1),key="-Nu2_fit-")],
    [sg.Text(text="Number of breakpoints:",font='Courier 18',size=(41,1)),
     sg.InputText("3",font='Courier 18',size=(3,1),key="-EN-",tooltip=">= 1"),
     sg.Text(text="Number of breakpoints:",font='Courier 18',size=(41,1)),
     sg.InputText("1",font='Courier 18',size=(3,1),key="-Nu1N-",tooltip=">= 0"),
     sg.Text(text="Number of breakpoints:",font='Courier 18',size=(39,1)),
     sg.InputText("1",font='Courier 18',size=(3,1),key="-Nu2N-",tooltip=">= 0")],
    [sg.Text("Data range:",font='Courier 18',size=(22,1),key="-ERangeText-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-E_start-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-E_stop-"),
     sg.Button("Re-Plot",font='Courier 18',key="-rePlotE-",size=(8,1)),
     sg.Text("Data range:",font='Courier 18',size=(23,1)),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-Nu1_start-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-Nu1_stop-"),
     sg.Button("Re-Plot",font='Courier 18',key="-rePlotNu1-",size=(7,1)),
     sg.Text("Data range:",font='Courier 18',size=(24,1)),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-Nu2_start-"),
     sg.InputText("0",font='Courier 18',size=(6,1),key="-Nu2_stop-"),
     sg.Button("Re-Plot",font='Courier 18',key="-rePlotNu2-",size=(7,1))],
    [sg.Button("Fit Elastic Modulus",font='Courier 18',size=(45,1),key="-fitElastic-",tooltip="Fitting may take several moments"),
     sg.Button("Fit first Poisson's ratio",font='Courier 18',size=(45,1),key="-fitNu1-",tooltip="Fitting may take several moments"),
     sg.Button("Fit second Poisson's ratio",font='Courier 18',size=(45,1),key="-fitNu2-",tooltip="Fitting may take several moments")],
    [sg.Image(key="-E_CANVAS-"),
     sg.Image(key="-Nu1_CANVAS-"),
     sg.Image(key="-Nu2_CANVAS-")],
    [sg.Text(text="E=\nYield strain = \nYield Strength = ",font='Courier 18',size=(45,3),key="-E_result-"),
     sg.Text("Nu1=",font='Courier 18',size=(45,3),key="-Nu1_result-"),
     sg.Text("Nu2=",font='Courier 18',size=(45,3),key="-Nu2_result-")],
    
]

# Create the Window
window = sg.Window("Uniaxial Analysis GUI", layout,finalize=True)

# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    print(event)
    if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks cancel
        break
    
    if event == '-LOAD-':
        # If you're loading a new file after completing an analysis,
        # clear the fitted figures
        if 'e_fig' in globals():
            plt.close(e_fig)
            window["-E_CANVAS-"].update()
        if 'Nu1_fig' in globals():
            plt.close(Nu1_fig)
            window["-Nu1_CANVAS-"].update()
        if 'Nu2_fig' in globals():
            plt.close(Nu2_fig)
            window["-Nu2_CANVAS-"].update()
        if 'fig_E' in globals():
            plt.close(fig_E)
            window["-E_CANVAS-"].update()
        if 'fig_Nu1' in globals():
            plt.close(fig_Nu1)
            window["-Nu1_CANVAS-"].update()
        if 'fig_Nu2' in globals():
            plt.close(fig_Nu2)
            window["-Nu2_CANVAS-"].update()
            
        log_file = values['-logfile-']
        if 'Comp' in log_file:
            compression_flag = True
        else:
            compression_flag = False
        os.chdir(values['-directory-'])
        strain_dir, primary_strain, primary_stress, \
            secondary_strain1, secondary_strain2, Nu1_label,\
            Nu2_label = uf.load_data(log_file)
        # Set parameter labels
        window['-E_fit-'].update(f"Fitting parameters for E{strain_dir}")
        window['-Nu1_fit-'].update(f"Fitting parameters for Nu{Nu1_label}")
        window['-Nu2_fit-'].update(f"Fitting parameters for Nu{Nu2_label}")
        window['-ERangeText-'].update(f"Max data points: {len(primary_strain)} ")
        # Plot to canvases
        fig_E, ax_E = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_E.scatter(primary_strain,primary_stress,color="grey",s=20)
        ax_E.set_xlabel(f'True strain in {strain_dir}')
        ax_E.set_ylabel(f'True stress in {strain_dir}, MPa')
        if compression_flag: # If a compressive uniaxial sim
            ax_E.invert_xaxis()
            # ax_E.invert_yaxis()
        draw_figure(window["-E_CANVAS-"], fig_E)
        
        fig_Nu1, ax_Nu1 = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_Nu1.scatter(primary_strain,secondary_strain1,color="grey",s=20)
        ax_Nu1.set_xlabel(f'True strain in {Nu1_label[0]}')
        ax_Nu1.set_ylabel(f'True strain in {Nu1_label[1]}')
        if compression_flag: # If a compressive uniaxial sim
            ax_Nu1.invert_xaxis()
            # ax_Nu1.invert_yaxis()
        draw_figure(window["-Nu1_CANVAS-"], fig_Nu1)
        
        fig_Nu2, ax_Nu2 = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_Nu2.scatter(primary_strain,secondary_strain2,color="grey",s=20)
        ax_Nu2.set_xlabel(f'True strain in {Nu2_label[0]}')
        ax_Nu2.set_ylabel(f'True strain in {Nu2_label[1]}')
        if compression_flag: # If a compressive uniaxial sim
            ax_Nu2.invert_xaxis()
            # ax_Nu2.invert_yaxis()
        draw_figure(window["-Nu2_CANVAS-"], fig_Nu2)
        window.refresh()
        
    if event == '-rePlotE-':
        ym_start = int(values["-E_start-"])
        ym_end = int(values["-E_stop-"])
        # Erase the non-fitted/fitted figure when replotting with the new range
        if 'e_fig' in globals():
            plt.close(e_fig)
            window["-E_CANVAS-"].update()
        else:
            plt.close(fig_E)
            window["-E_CANVAS-"].update()
            
        fig_E, ax_E = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_E.scatter(primary_strain[ym_start:ym_end],primary_stress[ym_start:ym_end],color="grey",s=20)
        ax_E.set_xlabel(f'True strain in {strain_dir}')
        ax_E.set_ylabel(f'True stress in {strain_dir}, GPa')
        if compression_flag:
            ax_E.invert_xaxis()
        draw_figure(window["-E_CANVAS-"], fig_E)
        window.refresh()
    
    if event == '-rePlotNu1-':
        # Erase the non-fitted/fitted figure when replotting with the new range
        if 'Nu1_fig' in globals():
            plt.close(Nu1_fig)
            window["-Nu1_CANVAS-"].update()
        else:
            plt.close(fig_Nu1)
            window["-Nu1_CANVAS-"].update()
            
        Nu1_start = int(values['-Nu1_start-'])
        Nu1_stop = int(values['-Nu1_stop-'])
        
        if Nu1_start < 0:
            Nu1_start = 0
            
        if Nu1_stop == 0:
            Nu1_stop = -1
            
        fig_Nu1, ax_Nu1 = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_Nu1.scatter(primary_strain[Nu1_start:Nu1_stop],secondary_strain1[Nu1_start:Nu1_stop],color="grey",s=20)
        ax_Nu1.set_xlabel(f'True strain in {Nu1_label[0]}')
        ax_Nu1.set_ylabel(f'True strain in {Nu1_label[1]}')
        if compression_flag:
            ax_Nu1.invert_xaxis()
        draw_figure(window["-Nu1_CANVAS-"], fig_Nu1)
        window.refresh()
        
    if event == '-rePlotNu2-':
        # Erase the non-fitted/fitted figure when replotting with the new range
        if 'Nu2_fig' in globals():
            plt.close(Nu2_fig)
            window["-Nu2_CANVAS-"].update()
        else:
            plt.close(fig_Nu2)
            window["-Nu2_CANVAS-"].update()
            
        Nu2_start = int(values['-Nu2_start-'])
        Nu2_stop = int(values['-Nu2_stop-'])
        
        if Nu2_start < 0:
            Nu2_start = 0
            
        if Nu2_stop == 0:
            Nu2_stop = -1
            
        fig_Nu2, ax_Nu2 = plt.subplots(figsize=(6.4,5),tight_layout=True)
        ax_Nu2.scatter(primary_strain[Nu2_start:Nu2_stop],secondary_strain2[Nu2_start:Nu2_stop],color="grey",s=20)
        ax_Nu2.set_xlabel(f'True strain in {Nu2_label[0]}')
        ax_Nu2.set_ylabel(f'True strain in {Nu2_label[1]}')
        if compression_flag:
            ax_Nu2.invert_xaxis()
        draw_figure(window["-Nu2_CANVAS-"], fig_Nu2)
        window.refresh()
        
    if event == '-fitElastic-' or event == '-FITALL-':
        # Erase the fitted figure when rerunning the analysis
        if 'e_fig' in globals():
            plt.close(e_fig)
            window["-E_CANVAS-"].update()
        else:
            plt.close(fig_E)
            window["-E_CANVAS-"].update()
            
        p_points = int(values["-EN-"])
        ym_start = int(values["-E_start-"])
        ym_end = int(values["-E_stop-"])
        primary_fit, youngsMod, ym_stderror, yield_strain, yield_strength\
            = uf.compute_elastic_modulus(primary_strain,primary_stress,p_points,data_range=(ym_start,ym_end),compression_flag=compression_flag)
        if 'sklearn' in str(type(primary_fit)):
            e_fig, e_ax = uf.plot_elastic(primary_fit, strain_dir, youngsMod, yield_strain, yield_strength,primary_strain,primary_stress)
            E_result_text = f"E{strain_dir} = {youngsMod:0.5f} GPa\nRsquared = {ym_stderror:0.6f}"
        else:
            e_fig, e_ax = uf.plot_elastic(primary_fit, strain_dir, youngsMod, yield_strain, yield_strength)
            E_result_text = f"E{strain_dir} = {youngsMod:0.5f} +- {ym_stderror:0.5f} GPa\nYield strain = {yield_strain:0.4f}\nYield strength = {yield_strength:0.2f} MPa"
        if compression_flag:
            e_ax.invert_xaxis()
        draw_figure(window["-E_CANVAS-"],e_fig)
        window["-E_result-"].update(E_result_text)
        
        window.refresh()
        
    if event == '-fitNu1-' or event == '-FITALL-':
        # Erase the fitted figure when rerunning the analysis
        if 'Nu1_fig' in globals():
            plt.close(Nu1_fig)
            window["-Nu1_CANVAS-"].update()
        else:
            plt.close(fig_Nu1)
            window["-Nu1_CANVAS-"].update()
            
        n1_points = int(values['-Nu1N-'])
        Nu1_start = int(values['-Nu1_start-'])
        Nu1_stop = int(values['-Nu1_stop-'])
        Nu1_fit, Nu1, Nu1_stderror\
            = uf.compute_poisson(primary_strain, secondary_strain1, n1_points,\
                                 Nu1_label,data_range=(Nu1_start,Nu1_stop))
        
        if 'sklearn' in str(type(Nu1_fit)):
            Nu1_fig, Nu1_ax = uf.plot_poisson(Nu1_fit, Nu1_label, Nu1,primary_strain,secondary_strain1)
            Nu1_result_text = f"Nu{Nu1_label} = {Nu1:0.6f}\nRsquared = {Nu1_stderror:0.6f}"
        else:
            Nu1_fig, Nu1_ax = uf.plot_poisson(Nu1_fit, Nu1_label, Nu1)
            Nu1_result_text = f"Nu{Nu1_label} = {Nu1:0.6f} +- {Nu1_stderror:0.6f}"
        if compression_flag:
            Nu1_ax.invert_xaxis()
        draw_figure(window["-Nu1_CANVAS-"],Nu1_fig)
        window["-Nu1_result-"].update(Nu1_result_text)
        window.refresh()
        
    if event == '-fitNu2-' or event == '-FITALL-':
        # Erase the fitted figure when rerunning the analysis
        if 'Nu2_fig' in globals():
            plt.close(Nu2_fig)
            window["-Nu2_CANVAS-"].update()
        else:
            plt.close(fig_Nu2)
            window["-Nu2_CANVAS-"].update()
            
        n2_points = int(values['-Nu2N-'])
        Nu2_start = int(values['-Nu2_start-'])
        Nu2_stop = int(values['-Nu2_stop-'])
        Nu2_fit, Nu2, Nu2_stderror\
            = uf.compute_poisson(primary_strain, secondary_strain2, n2_points,\
                                 Nu2_label,data_range=(Nu2_start,Nu2_stop))
        
        if 'sklearn' in str(type(Nu2_fit)):
            Nu2_fig, Nu2_ax = uf.plot_poisson(Nu2_fit, Nu2_label, Nu2,primary_strain,secondary_strain2)
            Nu2_result_text = f"Nu{Nu2_label} = {Nu2:0.6f}\nRsquared = {Nu2_stderror:0.6f}"
        else:
            Nu2_fig, Nu2_ax = uf.plot_poisson(Nu2_fit, Nu2_label, Nu2)
            Nu2_result_text = f"Nu{Nu2_label} = {Nu2:0.6f} +- {Nu2_stderror:0.6f}"
        
        if compression_flag:
            Nu2_ax.invert_xaxis()
        draw_figure(window["-Nu2_CANVAS-"],Nu2_fig)
        window["-Nu2_result-"].update(Nu2_result_text)
        window.refresh()

    if event == '-SAVE-' or event == '-FITALL-':
        uf.write_summaries(log_file,strain_dir,youngsMod,Nu1_label,Nu1,Nu2_label,\
                            Nu2,yield_strain,yield_strength,primary_fit,Nu1_fit,Nu2_fit,\
                            ym_stderror,Nu1_stderror,Nu2_stderror)
        primary_filename = log_file.split('.log')[0]+f'_stress_strain_{strain_dir}.pdf'
        Nu1_filename = log_file.split('.log')[0]+f'Nu{Nu1_label}.pdf'
        Nu2_filename = log_file.split('.log')[0]+f'Nu{Nu2_label}.pdf'
        e_fig.savefig(primary_filename,dpi=300)
        Nu1_fig.savefig(Nu1_filename,dpi=300)
        Nu2_fig.savefig(Nu2_filename,dpi=300)
        print(f"Data saved! Plots saved to {primary_filename}, {Nu1_filename}, and {Nu2_filename}.")

window.close()
