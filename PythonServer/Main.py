import os
from urllib.error import URLError

from PythonServer.Primer_Builder import *
from tkinter import *
import tkinter.messagebox


# ==================== On close ================================
def on_close_click():
    root.destroy()
    return


# ==================== Open output file on finish ================================

def open_output_file(specie_name, symbol_name):
    file = "cd ../output & notepad.exe primer_list_" + specie_name + "_" + symbol_name + ".txt"
    os.system(file)


# ==================== Back to default values ================================

def set_default_values():
    config_default = json.load(open("config defualt.json"))
    amplicon_max.set(config_default["Amplicon Length"]["Max"])
    amplicon_min.set(config_default["Amplicon Length"]["Min"])
    amplicon_opt.set(config_default["Amplicon Length"]["Avg"])
    leng_max.set(config_default["Length"]["Max"])
    leng_min.set(config_default["Length"]["Min"])
    leng_opt.set(config_default["Length"]["Avg"])
    tm_max.set(config_default["Tm"]["Max"])
    tm_min.set(config_default["Tm"]["Min"])
    tm_opt.set(config_default["Tm"]["Avg"])
    gc_max.set(config_default["GC Percent"]["Max"])
    gc_min.set(config_default["GC Percent"]["Min"])
    gc_opt.set(config_default["GC Percent"]["Avg"])
    tm_dif_max.set(config_default["Temperature difference"]["Max"])
    tm_dif_min.set(config_default["Temperature difference"]["Min"])


# ==================== Create finish view ================================

def finish_view(specie_name, symbol_name):
    root.geometry('{}x{}'.format(400, 100))
    center.destroy()
    specie_label_txt.set("Done! the primers file saved in the output folder.")
    symbol_label.destroy()
    entry_symbol.destroy()
    entry_specie.destroy()
    button_submit.destroy()
    button_close = Button(btm_frame, text="Close", width=20, command=on_close_click)
    button_file = Button(btm_frame, text="Open file", width=20,
                         command=lambda: open_output_file(specie_name, symbol_name))
    button_close.grid(row=0, column=0, padx=(40, 0))
    button_file.grid(row=0, column=1, padx=(5, 0))
    return


# ==================== error alert box ================================

def error_view(error):
    if error == "type":
        tkinter.messagebox.showerror("Error", "Wrong specie or symbol name")
    elif error == "url":
        tkinter.messagebox.showerror("Error", "No internet connection ")


# ==================== Submit, start primer building ================================

def on_submit_click():
    # Update config
    specie_name = specie.get()
    symbol_name = symbol.get()
    config["Amplicon Length"]["Max"] = amplicon_max.get()
    config["Amplicon Length"]["Min"] = amplicon_min.get()
    config["Amplicon Length"]["Avg"] = amplicon_opt.get()
    config["Length"]["Max"] = leng_max.get()
    config["Length"]["Min"] = leng_min.get()
    config["Length"]["Avg"] = leng_opt.get()
    config["Tm"]["Max"] = tm_max.get()
    config["Tm"]["Min"] = tm_min.get()
    config["Tm"]["Avg"] = tm_opt.get()
    config["GC Percent"]["Max"] = gc_max.get()
    config["GC Percent"]["Min"] = gc_min.get()
    config["GC Percent"]["Avg"] = gc_min.get()
    config["Temperature difference"]["Max"] = tm_dif_max.get()
    config["Temperature difference"]["Min"] = tm_dif_min.get()
    jsonFile = open("config.json", "w+")
    jsonFile.write(json.dumps(config))
    jsonFile.close()

    # Start primer builder
    try:
        build(specie_name, symbol_name)
        finish_view(specie_name, symbol_name)
    except TypeError:
        error_view("type")
    except URLError:
        error_view("url")
    return


# ==================== main , init GUI user interface ================================

if __name__ == '__main__':
    config = json.load(open("config.json"))
    root = Tk()
    root.title('Primer Builder')
    root.geometry('{}x{}'.format(570, 300))
    # create all of the main containers
    top_frame = Frame(root, bg='red', width=450, height=50, pady=3)
    center = Frame(root, width=490, height=20, padx=3, pady=3)
    btm_frame = Frame(root, bg='white', width=450, height=45, pady=3)

    # layout all of the main containers
    root.grid_rowconfigure(1, weight=1)
    root.grid_columnconfigure(0, weight=1)

    top_frame.grid(row=0, sticky="ew")
    center.grid(row=1, sticky="nsew")
    btm_frame.grid(row=3, sticky="ew")

    # create the widgets for the top frame
    main_label = StringVar()
    main_label.set("Primer Builder")
    specie = StringVar()
    specie.set("")
    symbol = StringVar()
    symbol.set("")
    specie_label_txt = StringVar()
    specie_label_txt.set("specie:")
    primer_label = Label(top_frame, textvariable=main_label, bg="yellow")
    specie_label = Label(top_frame, textvariable=specie_label_txt)
    symbol_label = Label(top_frame, text='Symbol:')
    entry_specie = Entry(top_frame, textvariable=specie)
    entry_symbol = Entry(top_frame, textvariable=symbol)

    # layout the widgets in the top frame
    primer_label.grid(row=0, columnspan=100, padx=(80, 0), pady=10)
    specie_label.grid(row=1, column=0, padx=(80, 10))
    symbol_label.grid(row=1, column=2, padx=10)
    entry_specie.grid(row=1, column=1, padx=2)
    entry_symbol.grid(row=1, column=3, padx=2)

    amplicon_max = StringVar()
    amplicon_max.set(config["Amplicon Length"]["Max"])
    amplicon_min = StringVar()
    amplicon_min.set(config["Amplicon Length"]["Min"])
    amplicon_opt = StringVar()
    amplicon_opt.set(config["Amplicon Length"]["Avg"])
    amplicon_max_label = Label(center, text='Amplicon max length:')
    amplicon_min_label = Label(center, text='Amplicon min length:')
    amplicon_opt_label = Label(center, text='Amplicon optimal length:')
    entry_amp_max = Entry(center, textvariable=amplicon_max)
    entry_amp_min = Entry(center, textvariable=amplicon_min)
    entry_amp_opt = Entry(center, textvariable=amplicon_opt)

    leng_max = StringVar()
    leng_max.set(config["Length"]["Max"])
    leng_min = StringVar()
    leng_min.set(config["Length"]["Min"])
    leng_opt = StringVar()
    leng_opt.set(config["Length"]["Avg"])
    leng_max_label = Label(center, text='Primer max length:')
    leng_min_label = Label(center, text='Primer min length:')
    leng_opt_label = Label(center, text='Primer optimal length:')
    entry_len_max = Entry(center, textvariable=leng_max)
    entry_len_min = Entry(center, textvariable=leng_min)
    entry_len_opt = Entry(center, textvariable=leng_opt)

    tm_max = StringVar()
    tm_max.set(config["Tm"]["Max"])
    tm_min = StringVar()
    tm_min.set(config["Tm"]["Min"])
    tm_opt = StringVar()
    tm_opt.set(config["Tm"]["Avg"])
    tm_max_label = Label(center, text='Primer max Tm:')
    tm_min_label = Label(center, text='Primer min Tm:')
    tm_opt_label = Label(center, text='Primer optimal Tm:')
    entry_tm_max = Entry(center, textvariable=tm_max)
    entry_tm_min = Entry(center, textvariable=tm_min)
    entry_tm_opt = Entry(center, textvariable=tm_opt)

    gc_max = StringVar()
    gc_max.set(config["GC Percent"]["Max"])
    gc_min = StringVar()
    gc_min.set(config["GC Percent"]["Min"])
    gc_opt = StringVar()
    gc_opt.set(config["GC Percent"]["Avg"])
    gc_max_label = Label(center, text='Primer max GC%:')
    gc_min_label = Label(center, text='Primer min GC%:')
    gc_opt_label = Label(center, text='Primer optimal GC%:')
    entry_gc_max = Entry(center, textvariable=gc_max)
    entry_gc_min = Entry(center, textvariable=gc_min)
    entry_gc_opt = Entry(center, textvariable=gc_opt)

    tm_dif_max = StringVar()
    tm_dif_max.set(config["Temperature difference"]["Max"])
    tm_dif_min = StringVar()
    tm_dif_min.set(config["Temperature difference"]["Min"])
    tm_dif_max_label = Label(center, text='Primers Tm max dif:')
    tm_dif_min_label = Label(center, text='Primer Tm min dif:')
    entry_tm_dif_max = Entry(center, textvariable=tm_dif_max)
    entry_tm_dif_min = Entry(center, textvariable=tm_dif_min)

    amplicon_min_label.grid(row=0, column=0, padx=5)
    amplicon_max_label.grid(row=1, column=0, padx=5)
    amplicon_opt_label.grid(row=2, column=0, padx=5)
    leng_min_label.grid(row=3, column=0, padx=5)
    leng_max_label.grid(row=4, column=0, padx=5)
    leng_opt_label.grid(row=5, column=0, padx=5)
    tm_min_label.grid(row=6, column=0, padx=5)
    tm_max_label.grid(row=7, column=0, padx=5)
    tm_opt_label.grid(row=8, column=0, padx=5)
    gc_min_label.grid(row=0, column=2, padx=5)
    gc_max_label.grid(row=1, column=2, padx=5)
    gc_opt_label.grid(row=2, column=2, padx=5)
    tm_dif_min_label.grid(row=3, column=2, padx=5)
    tm_dif_max_label.grid(row=4, column=2, padx=5)

    entry_amp_min.grid(row=0, column=1)
    entry_amp_max.grid(row=1, column=1)
    entry_amp_opt.grid(row=2, column=1)
    entry_len_min.grid(row=3, column=1)
    entry_len_max.grid(row=4, column=1)
    entry_len_opt.grid(row=5, column=1)
    entry_tm_min.grid(row=6, column=1)
    entry_tm_max.grid(row=7, column=1)
    entry_tm_opt.grid(row=8, column=1)
    entry_gc_min.grid(row=0, column=3)
    entry_gc_max.grid(row=1, column=3)
    entry_gc_opt.grid(row=2, column=3)
    entry_tm_dif_min.grid(row=3, column=3)
    entry_tm_dif_max.grid(row=4, column=3)

    button_default = Button(center, text="Default values", width=20, command=set_default_values)
    button_default.grid(row=5, column=2)

    button_submit = Button(btm_frame, text="Get Primers", width=20, command=on_submit_click)
    button_default = Button(btm_frame, text="Default values", width=20, command=set_default_values)
    button_submit.grid(row=0, column=0, padx=(195, 0))
    root.mainloop()