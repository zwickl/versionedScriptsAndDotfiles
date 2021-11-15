import sys
from Tkinter import *
import tkFileDialog
import tkFont
from ttk import *
import argparse
from textwrap import fill
import re
import shlex
from plotutils import ArgparseActionAppendToDefault

widget_padx = 10
text_widths = 60

class ArgparseOption(object):
    def __init__(self, option):
        #use the last listed flag, which is likely to be the more descriptive long one
        #there will be no option_strings for a positional arg
        if option.option_strings:
            self.output_arg = option.option_strings[-1]
        else:
            self.output_arg = None
        
        if 'HIDE' in option.help:
            self.hide = True
            self.return_string = []
        else:
            self.hide = False
            self.return_string = []

        self.nargs = option.nargs


class ArgparseBoolOption(ArgparseOption):
    def __init__(self, option, frame, row=-1, column=0):
        ArgparseOption.__init__(self, option)
       
        #using Variable rather than IntVar since it allows a default of None
        self.var = Variable()
        self.var.set(option.default)

        if not self.hide:
            help_string = re.sub('[(]default [)]', '', option.help).strip()

            if help_string:
                label_string = help_string
            else:
                label_string = re.sub('--', '', option.option_strings[-1])

            Label(frame, text=fill(label_string, text_widths)).grid(row=row, column=column)
            
            self.widget = Checkbutton(frame, variable=self.var, onvalue=1, offvalue=0)
            '''
            if isinstance(option, argparse._StoreTrueAction):
                self.var.set(False)
            else:
                self.var.set(True)
            '''
            
            if isinstance(option, argparse._StoreTrueAction):
                self.widget = Checkbutton(frame, variable=self.var, onvalue=1, offvalue=0)
            elif isinstance(option, argparse._StoreFalseAction):
                self.widget = Checkbutton(frame, variable=self.var, onvalue=0, offvalue=1)
            
            if row < 0:
                self.widget.pack()
            else:
                self.widget.grid(row=row, column=column+1, padx=widget_padx)

    def make_string(self):
        #print self
        if bool(self.var.get()):
            return [ self.output_arg ]
        else:
            return []


class ArgparseStringOption(ArgparseOption):
    def __init__(self, option, frame, row=-1, column=0):
        ArgparseOption.__init__(self, option)
        self.var = StringVar()
        self.widget = Entry(frame, textvariable=self.var)
        if option.required:
            req_string = 'REQ: '
        else:
            req_string = ''

        if not self.hide:
            help_string = re.sub('[(]default [)]', '', option.help).strip()
            if help_string:
                label_string = req_string + help_string
            else:
                label_string = req_string + re.sub('--', '', option.option_strings[-1])

            if row < 0:
                Label(frame, text=fill(label_string, text_widths)).pack(side=LEFT)
                self.widget.pack()
            else:
                Label(frame, text=fill(label_string, text_widths)).grid(row=row, column=column, padx=widget_padx)
                self.widget.grid(row=row, column=column+1, padx=widget_padx)
        
        if option.default:
            if isinstance(option.default, list):
                self.widget.insert(0, ' '.join([str(val) for val in option.default]))
            else:
                self.widget.insert(0, str(option.default))


    def make_string(self):
        #print self
        if self.var.get():
            #print self.output_arg, self.var.get(), type(self.var.get()), self.nargs
            self.return_string.append(self.output_arg)
            
            if self.nargs and (self.nargs in [ '*', '+' ] or self.nargs > 1):
                #shlex.split here properly leaves quoted strings unsplit
                splt = shlex.split(self.var.get())
                #this is an annoying special case, where a leading "-" in an argument has to 
                #have double quotes explicitly embedded in the string
                for num, s in enumerate(splt):
                    if s[0] == '-':
                        splt[num] = '"' + s + '"'

                self.return_string.extend(splt)
            else:
                self.return_string.append(self.var.get())
        #print '\t', return_string 
        return self.return_string


class ArgparseOptionMenuOption(ArgparseOption):
    def __init__(self, option, frame, row=-1, column=0):
        #print option
        self.var = StringVar()

        if option.required:
            req_string = 'REQ: '
        else:
            req_string = ''
        
        self.var.set(option.default)

        #OptionMenu signature is this:
        #__init__(self, master, variable, value, *values, **kwargs)
        #where variable is "the resource textvariable", and value is the 
        #default value
        self.widget = OptionMenu(frame, self.var, option.choices[0], *option.choices)

        help_string = re.sub('[(]default [)]', '', option.help).strip()
        if help_string:
            label_string = req_string + help_string
        else:
            label_string = req_string + re.sub('--', '', option.option_strings[-1])

        if row < 0:
            Label(frame, text=fill(label_string, text_widths)).pack(side=LEFT)
            self.widget.pack()
        else:
            Label(frame, text=fill(label_string, text_widths)).grid(row=row, column=column, padx=widget_padx)
            self.widget.grid(row=row, column=column+1, padx=widget_padx)
 
        #Label(frame, text=fill(req_string + option.help, text_widths)).grid(row=row, column=column, padx=widget_padx)
        #self.widget.grid(row=row, column=column+1, padx=widget_padx)
       
        ArgparseOption.__init__(self, option)

    def make_string(self):
        #print self
        if self.var.get():
            self.return_string.append(self.output_arg)
            self.return_string.append(self.var.get())
        return self.return_string


class ArgparseFileOption(ArgparseOption):
    def __init__(self, option, frame, row=-1, column=0):
        ArgparseOption.__init__(self, option)

        Label(frame, text=fill(option.help)).grid(row=row, column=column, padx=widget_padx)
        if option.type and hasattr(option.type, "_mode"):
            #a mode would be here if the option is specified a file to argparse, rather than the path to a file
            if 'r' in option.type._mode:
                if self.nargs and (self.nargs in [ '*', '+' ] or self.nargs > 1):
                    Button(frame, text='OPEN', command=self.open_multiple_files_dialog).grid(row=row, column=column+1, padx=widget_padx) 
                else:
                    Button(frame, text='OPEN', command=self.open_file_dialog).grid(row=row, column=column+1, padx=widget_padx) 
            elif 'w' in option.type._mode:
                Button(frame, text='OPEN', command=self.output_file_dialog).grid(row=row, column=column+1, padx=widget_padx) 
        else:
            #this is obviously a total hack, and depends on the "destination" variable name assigned in argparse
            if 'out' in option.dest.lower():
                Button(frame, text='SAVE AS', command=self.output_file_dialog).grid(row=row, column=column+1, padx=widget_padx) 
            else:
                if self.nargs and (self.nargs in [ '*', '+' ] or self.nargs > 1):
                    Button(frame, text='OPEN', command=self.open_multiple_files_dialog).grid(row=row, column=column+1, padx=widget_padx) 
                else:
                    Button(frame, text='OPEN', command=self.open_file_dialog).grid(row=row, column=column+1, padx=widget_padx) 

        self.var = None

    def open_file_dialog(self):
        self.var = tkFileDialog.askopenfilename()

    def open_multiple_files_dialog(self):
        self.var = tkFileDialog.askopenfilenames()

    def output_file_dialog(self):
        self.var = tkFileDialog.asksaveasfilename()

    def make_string(self):
        #print self
        if self.var:
            if self.output_arg:
                self.return_string.append(self.output_arg)
            #self.var is actually a tuple here in the case of multiple filenames,
            #but may be just a string otherwise
            if isinstance(self.var, tuple):
                self.var = list(self.var)
            if isinstance(self.var, list):
                self.return_string.extend(self.var)
            else:
                self.return_string.append(self.var)
        return self.return_string


class ArgparseGui(object):
    #def __init__(self, tk, parser, height=768, width=1024):
    def __init__(self, parser, tk=None, height=768, width=1024):
        self.tk = tk or Tk()

        self.max_widgets_per_column = 20
        self.column_offset = 0
        self.row = 0
        
        ###############
        #from http://stackoverflow.com/questions/3085696/adding-a-scrollbar-to-a-grid-of-widgets-in-tkinter
        self.canvas = Canvas(self.tk, height=height, width=width, borderwidth=0, background="#ffffff")
        
        #background used to work here. wtf? 
        #self.frame = Frame(self.canvas, background="#ffffff")
        self.frame = Frame(self.canvas)
        
        self.vsb = Scrollbar(self.tk, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.pack(side="right", fill="y")
        #self.canvas.pack(side="left", fill="both", expand=True)
        #self.canvas.create_window((4,4), window=self.frame, anchor="nw", tags="self.frame")
        #self.frame.bind("<Configure>", self.OnFrameConfigure)


        self.hsb = Scrollbar(self.tk, orient="horizontal", command=self.canvas.xview)
        self.canvas.configure(xscrollcommand=self.hsb.set)

        self.hsb.pack(side="bottom", fill="x")
        
        
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((4,4), window=self.frame, anchor="nw", 
                                      tags="self.frame")
        self.frame.bind("<Configure>", self.OnFrameConfigure)
        ################

        #self.frame = Frame(tk)
        #self.frame.pack(expand=False)
        
        self.option_list = []
        #bizarrely options appear once for each potential flag (i.e. short and long), so need to keep track
        seenOptions = []

        '''
        group_list = [parser._action_groups[0]]
        if len(parser._action_groups) > 2:
            group_list.extend(parser._action_groups[2:])
        group_list.append(parser._action_groups[1])
        for group in group_list:
            if group.title != "optional arguments":
                display_title = group.title
            else:
                display_title = "Misc. options" 
            Label(self.frame, text=fill(display_title)).grid(row=5)
            for option in group._group_actions:
                if option not in seenOptions:
                    seenOptions.append(option)
                    if isinstance(option, argparse._StoreTrueAction):
                        gui_option = ArgparseBoolOption(option, self.frame, len(self.option_list))
                        self.option_list.append(gui_option)
                    elif isinstance(option, argparse._StoreAction):
                        if isinstance(option.type, type(str)):
                            gui_option = ArgparseStringOption(option, self.frame, len(self.option_list))
                            self.option_list.append(gui_option)
                        elif isinstance(option.type, argparse.FileType):
                            gui_option = ArgparseFileOption(option, self.frame, len(self.option_list))
                            self.option_list.append(gui_option)
        '''
        
        #for optName, option in parser._option_string_actions.items():
        group_list = [parser._action_groups[0]]
        if len(parser._action_groups) > 2:
            group_list.extend(parser._action_groups[2:])
        group_list.append(parser._action_groups[1])
        for group in group_list:
            #print
            #print group.title
            if len(group._group_actions):
                group_title_displayed = False
                for option in group._group_actions:
                    if option not in seenOptions:
                        #print '\t', option
                        if not group_title_displayed:
                            if self.row + len(group._group_actions) > self.max_widgets_per_column and self.row > self.max_widgets_per_column - 3:
                                self.row = 0
                                self.column_offset += 2
                            if group.title != "optional arguments":
                                display_title = group.title.upper()
                            else:
                                display_title = "Misc. options".upper()
                            Label(self.frame, text=fill(display_title), font=tkFont.Font(size=14, weight='bold')).grid(row=self.row, column=self.column_offset, columnspan=2)
                            group_title_displayed = True
                            self.row += 1
                        
                        seenOptions.append(option)
                        if isinstance(option, argparse._StoreTrueAction) or isinstance(option, argparse._StoreFalseAction):
                            gui_option = ArgparseBoolOption(option, self.frame, row=self.row, column=self.column_offset)
                        elif isinstance(option, argparse._StoreAction) or isinstance(option, argparse._AppendAction):
                            if option.choices:
                                gui_option = ArgparseOptionMenuOption(option, self.frame, row=self.row, column=self.column_offset)
                            elif (isinstance(option.type, type(str)) or option.type is None) and 'file' in option.dest.lower():
                                gui_option = ArgparseFileOption(option, self.frame, row=self.row, column=self.column_offset)
                            #if no type is specified to ArgumentParser.add_argument then the default is str
                            elif option.type is None or isinstance(option.type, type(str)):
                                gui_option = ArgparseStringOption(option, self.frame, row=self.row, column=self.column_offset)
                            elif isinstance(option.type, argparse.FileType):
                                gui_option = ArgparseFileOption(option, self.frame, row=self.row, column=self.column_offset)
                            else:
                                print "unknown Store action:", option
                                sys.exit()
                        elif isinstance(option, ArgparseActionAppendToDefault):
                            gui_option = ArgparseStringOption(option, self.frame, row=self.row, column=self.column_offset)
                        elif isinstance(option, argparse._HelpAction):
                            continue
                        else:
                            print "unknown action:", option
                            sys.exit()
                        
                        self.row += 1
                        if self.row >= self.max_widgets_per_column:
                            self.row = 0
                            self.column_offset += 2

                        self.option_list.append(gui_option)

        done = Button(self.frame, text='DONE', command=self.done).grid(row=self.max_widgets_per_column+1, column=0, padx=widget_padx) 
        #done.focus_set()
        cancel = Button(self.frame, text='Cancel', command=self.cancel).grid(row=self.max_widgets_per_column+1, column=1, padx=widget_padx) 
        #Button(self.frame, text='DONE', command=self.done).grid(row=self.row, column=self.column_offset, padx=widget_padx) 
        #Button(self.frame, text='Cancel', command=self.cancel).grid(row=self.row, column=self.column_offset + 1, padx=widget_padx) 

        self.cancelled = False


    def make_commandline_list(self):
        #print self
        return_list = []
        for option in self.option_list:
            return_list.extend(option.make_string())
        print return_list
        return return_list

    def done(self):
        self.frame.destroy()

    def cancel(self):
        self.frame.destroy()
        self.cancelled = True

    def OnFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

