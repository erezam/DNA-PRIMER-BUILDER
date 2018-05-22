
from PythonServer.Primer_Builder import *
from Tkinter import *


def on_close_click():
    app.destroy()
    return

def set_default_values():
    config_default = json.load(open("config defualt.json"))
    amplicon_max.set(config_default["Amplicon Length"]["Max"])
    amplicon_min.set(config_default["Amplicon Length"]["Min"])
    leng_max.set(config_default["Length"]["Max"])
    leng_min.set(config_default["Length"]["Min"])
    tm_max.set(config_default["Tm"]["Max"])
    tm_min.set(config_default["Tm"]["Min"])
    gc_max.set(config_default["GC Percent"]["Max"])
    gc_min.set(config_default["GC Percent"]["Min"])
    tm_dif_max.set(config_default["Temperature difference"]["Max"])
    tm_dif_min.set(config_default["Temperature difference"]["Min"])


def on_submit_click():
    # Update config
    config["Amplicon Length"]["Max"] = amplicon_max.get()
    config["Amplicon Length"]["Min"] = amplicon_min.get()
    config["Length"]["Max"] = leng_max.get()
    config["Length"]["Min"] = leng_min.get()
    config["Tm"]["Max"] = tm_max.get()
    config["Tm"]["Min"] = tm_min.get()
    config["GC Percent"]["Max"] = gc_max.get()
    config["GC Percent"]["Min"] = gc_min.get()
    config["Temperature difference"]["Max"] = tm_dif_max.get()
    config["Temperature difference"]["Min"] = tm_dif_min.get()
    jsonFile = open("config.json", "w+")
    jsonFile.write(json.dumps(config))
    jsonFile.close()

    # Start primer builder
    transcript_data(specie.get(), symbol.get())
    return


if __name__ == '__main__':
    config = json.load(open("config.json"))
    root = Tk()
    root.title('Primer Builder')
    root.geometry('{}x{}'.format(500, 250))
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
    primer_label = Label(top_frame, textvariable=main_label, bg="yellow")
    specie_label = Label(top_frame, text='Specie:')
    symbol_label = Label(top_frame, text='Symbol:')
    entry_specie = Entry(top_frame, textvariable=specie)
    entry_symbol = Entry(top_frame, textvariable=symbol)

    # layout the widgets in the top frame
    primer_label.grid(row=0, columnspan=100, padx=(50, 0), pady=10)
    specie_label.grid(row=1, column=0, padx=(50, 10))
    symbol_label.grid(row=1, column=2, padx=10)
    entry_specie.grid(row=1, column=1, padx=2)
    entry_symbol.grid(row=1, column=3, padx=2)

    amplicon_max = StringVar()
    amplicon_max.set(config["Amplicon Length"]["Max"])
    amplicon_min = StringVar()
    amplicon_min.set(config["Amplicon Length"]["Min"])
    amplicon_max_label = Label(center, text='Amplicon max length:')
    amplicon_min_label = Label(center, text='Amplicon min length:')
    entry_amp_max = Entry(center, textvariable=amplicon_max)
    entry_amp_min = Entry(center, textvariable=amplicon_min)

    leng_max = StringVar()
    leng_max.set(config["Length"]["Max"])
    leng_min = StringVar()
    leng_min.set(config["Length"]["Min"])
    leng_max_label = Label(center, text='Primer max length:')
    leng_min_label = Label(center, text='Primer min length:')
    entry_len_max = Entry(center, textvariable=leng_max)
    entry_len_min = Entry(center, textvariable=leng_min)

    tm_max = StringVar()
    tm_max.set(config["Tm"]["Max"])
    tm_min = StringVar()
    tm_min.set(config["Tm"]["Min"])
    tm_max_label = Label(center, text='Primer max Tm:')
    tm_min_label = Label(center, text='Primer min Tm:')
    entry_tm_max = Entry(center, textvariable=tm_max)
    entry_tm_min = Entry(center, textvariable=tm_min)

    gc_max = StringVar()
    gc_max.set(config["GC Percent"]["Max"])
    gc_min = StringVar()
    gc_min.set(config["GC Percent"]["Min"])
    gc_max_label = Label(center, text='Primer max GC%:')
    gc_min_label = Label(center, text='Primer min GC%:')
    entry_gc_max = Entry(center, textvariable=gc_max)
    entry_gc_min = Entry(center, textvariable=gc_min)

    tm_dif_max = StringVar()
    tm_dif_max.set(config["Temperature difference"]["Max"])
    tm_dif_min = StringVar()
    tm_dif_min.set(config["Temperature difference"]["Min"])
    tm_dif_max_label = Label(center, text='Primers Tm max dif:')
    tm_dif_min_label = Label(center, text='Primer Tm max dif:')
    entry_tm_dif_max = Entry(center, textvariable=tm_dif_max)
    entry_tm_dif_min = Entry(center, textvariable=tm_dif_min)

    amplicon_min_label.grid(row=0, column=0, padx=5)
    amplicon_max_label.grid(row=1, column=0, padx=5)
    leng_min_label.grid(row=2, column=0, padx=5)
    leng_max_label.grid(row=3, column=0, padx=5)
    tm_min_label.grid(row=4, column=0, padx=5)
    tm_max_label.grid(row=5, column=0, padx=5)
    gc_min_label.grid(row=0, column=2, padx=5)
    gc_max_label.grid(row=1, column=2, padx=5)
    tm_dif_min_label.grid(row=2, column=2, padx=5)
    tm_dif_max_label.grid(row=3, column=2, padx=5)

    entry_amp_min.grid(row=0, column=1)
    entry_amp_max.grid(row=1, column=1)
    entry_len_min.grid(row=2, column=1)
    entry_len_max.grid(row=3, column=1)
    entry_tm_min.grid(row=4, column=1)
    entry_tm_max.grid(row=5, column=1)
    entry_gc_min.grid(row=0, column=3)
    entry_gc_max.grid(row=1, column=3)
    entry_tm_dif_min.grid(row=2, column=3)
    entry_tm_dif_max.grid(row=3, column=3)

    button_default = Button(center, text="Default values", width=20, command=set_default_values)
    button_default.grid(row=5, column=2)

    button_submit = Button(btm_frame, text="Get Primers", width=20, command=on_submit_click)
    button_default = Button(btm_frame, text="Default values", width=20, command=set_default_values)
    button_submit.grid(row=0, column=0, padx=(170, 0))
    root.mainloop()