
from PythonServer.Primer_Builder import *
from Tkinter import *


def on_close_click():
    app.destroy()
    return


def on_submit_click():
    transcript_data(specie.get(), symbol.get())
    name = "Done! Proceed to the output folder to view the results."
    label_specie_text.set(name)
    specieInput.destroy()
    label_symbol.destroy()
    symbolInput.destroy()
    button.configure(text="Ok", command=on_close_click)
    return


if __name__ == '__main__':
    root = Tk()
    root.title('Model Definition')
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
    primer_label = Label(top_frame, text='Primer Builder', bg="yellow")
    specie_label = Label(top_frame, text='Specie:')
    symbol_label = Label(top_frame, text='Symbol:')
    entry_specie = Entry(top_frame)
    entry_symbol = Entry(top_frame)

    # layout the widgets in the top frame
    primer_label.grid(row=0, columnspan=100, padx=(50,0), pady=10)
    specie_label.grid(row=1, column=0, padx=(50, 10))
    symbol_label.grid(row=1, column=2, padx=10)
    entry_specie.grid(row=1, column=1, padx=2)
    entry_symbol.grid(row=1, column=3, padx=2)

    # create the center widgets
    #center.grid_rowconfigure(0, weight=1)
    #center.grid_columnconfigure(1, weight=1)

    amplicon_max_label = Label(center, text='Amplicon max length:')
    amplicon_min_label = Label(center, text='Amplicon min length:')
    entry_amp_max = Entry(center)
    entry_amp_min = Entry(center)

    leng_max_label = Label(center, text='Primer max length:')
    leng_min_label = Label(center, text='Primer min length:')
    entry_len_max = Entry(center)
    entry_len_min = Entry(center)

    tm_max_label = Label(center, text='Primer max Tm:')
    tm_min_label = Label(center, text='Primer min Tm:')
    entry_tm_max = Entry(center)
    entry_tm_min = Entry(center)

    gc_max_label = Label(center, text='Primer max GC%:')
    gc_min_label = Label(center, text='Primer min GC%:')
    entry_gc_max = Entry(center)
    entry_gc_min = Entry(center)

    tm_dif_max_label = Label(center, text='Primer max Tm:')
    tm_dif_min_label = Label(center, text='Primer min Tm:')
    entry_tm_dif_max = Entry(center)
    entry_tm_dif_min = Entry(center)

    amplicon_max_label.grid(row=0, column=0, padx=5)
    amplicon_min_label.grid(row=1, column=0, padx=5)
    leng_max_label.grid(row=2, column=0, padx=5)
    leng_min_label.grid(row=3, column=0, padx=5)
    tm_max_label.grid(row=4, column=0, padx=5)
    tm_min_label.grid(row=5, column=0, padx=5)
    gc_max_label.grid(row=0, column=2, padx=5)
    gc_min_label.grid(row=1, column=2, padx=5)
    tm_dif_max_label.grid(row=2, column=2, padx=5)
    tm_dif_min_label.grid(row=3, column=2, padx=5)

    entry_amp_max.grid(row=0, column=1)
    entry_amp_min.grid(row=1, column=1)
    entry_len_max.grid(row=2, column=1)
    entry_len_min.grid(row=3, column=1)
    entry_tm_max.grid(row=4, column=1)
    entry_tm_min.grid(row=5, column=1)
    entry_gc_max.grid(row=0, column=3)
    entry_gc_min.grid(row=1, column=3)
    entry_tm_dif_max.grid(row=2, column=3)
    entry_tm_dif_min.grid(row=3, column=3)

    button = Button(btm_frame, text="Get Primers", width=20, command=on_submit_click)
    button.grid(row=0, column=0, padx=(170, 0))
    root.mainloop()

    '''
    app = Tk()
    app.title("Primer Builder")

    label_specie_text = StringVar()
    label_specie_text.set("Enter specie:")
    label_specie = Label(app, textvariable=label_specie_text, height=1)
    label_specie.pack(side=LEFT)

    specie = StringVar(None)
    specieInput = Entry(app, textvariable=specie)
    specieInput.pack(side=LEFT, pady=5)

    label_symbol_text = StringVar()
    label_symbol_text.set("Enter gene symbol")
    label_symbol = Label(app, textvariable=label_symbol_text, height=1)
    label_symbol.pack()

    symbol = StringVar(None)
    symbolInput = Entry(app, textvariable=symbol)
    symbolInput.pack(side=LEFT,pady=5)

    button = Button(app, text="Get Primers", width=20, command=on_submit_click)
    button.pack(side='bottom', padx=10, pady=10)

    # Config settings
    label_amplicon_max = Label(app, text="Amplicon max length:", height=1, width=50)
    label_amplicon_max.pack()
    amplicon_max = StringVar(None)
    amplicon_max_input = Entry(app, textvariable=amplicon_max)
    amplicon_max_input.pack(pady=5)


    app.mainloop()'''