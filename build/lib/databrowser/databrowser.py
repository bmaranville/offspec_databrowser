#!/usr/bin/python
import sys, os
#sys.path.append(os.path.join(os.environ['HOME'],'bin'))
import wx
import wx.grid
import wx.lib.masked
import wx.calendar
#from IPython.gui.wx import ipython_view
import reflectometry.reduction as red
from glob import glob1
#import cm
#from angle_bmp import *
#from plot_2d3 import *
from matplotlib.image import FigureImage
from matplotlib.figure import Figure
from matplotlib.pyplot import figure
from matplotlib.axis import XAxis, YAxis
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib import cm
from numpy import log
import numpy
import hashlib
import time, datetime
from .wxrebin_polcorr import rebinned_data, Supervisor
from .he3analyzer import He3Analyzer 

import xml.dom.minidom
import simplejson

class DataBrowser(wx.Frame):
    """ DataBrowser tool.  Loads all binned 2d data in directory, and allows joining of files
    Generates a file_catalog.txt file with information on created data objects """
    def __init__(self, data_directory = None, parent = None, load_saved_objects = True):
        wx.Frame.__init__(self, parent, title="MAGIK File Browser", size=(1024,768))
        #self.file_objects = file_objects
        self.data_directory = data_directory
        if not self.data_directory:
            self.data_directory = os.getcwd()
            print 'directory: ', os.getcwd()
        self.file_catalog = xml.dom.minidom.Document()
        self.file_catalog_path = os.path.join(self.data_directory, 'file_catalog.xml')
        file_objects, file_list, user_descriptions, elapsed_time = self.ProcessDirectory()
        self.SaveFileCatalog()
        self.file_objects = file_objects
        self.file_list = file_list
        self.obj_catalog = xml.dom.minidom.Document()
        obj_catalog_docroot = self.obj_catalog.createElement('object_catalog')
        self.obj_catalog.appendChild(obj_catalog_docroot)
        self.obj_catalog_path = os.path.join(self.data_directory, 'obj_catalog.xml')
        self.He3_cells = {} # start with no cells defined
        
        # replicate Supervisor functionality
        self.rb_count = 0
        self.rebinned_data_objects = []
        self.data_objects = {}
        self.plottable_count = 0
        self.plottable_2d_data_objects = []
        self.rebinned_data_treeitem = {}
        self.plottable_2d_data_treeitem = {}

        
        print "initializing frame"
        menubar = wx.MenuBar()
        filemenu = wx.Menu()
        quit = filemenu.Append(wx.ID_EXIT, '&Quit\tCtrl+Q', "Close window")
        self.Bind(wx.EVT_MENU, self.onExit, quit)
        #quit.SetBitmap(wx.Bitmap('icons/exit.png'))
        #filemenu.AppendItem(quit)
        save_csv = filemenu.Append(wx.ID_ANY, 'Save grid to &CSV (spreadsheet) file', 'Save grid to spreadsheet file')
        self.Bind(wx.EVT_MENU, self.SaveCatalogSpreadsheet, save_csv)
        imp_ts = filemenu.Append(wx.ID_ANY, '&Get elapsed time from original file timestamp', 'Get elapsed time from original file timestamp')
        self.Bind(wx.EVT_MENU, self.GetTimestamps, imp_ts)
        
        self.cellmenu = wx.Menu()
        add_cell = self.cellmenu.Append(wx.ID_ANY, '&Add new cell to list', 'Add new cell to list')
        self.Bind(wx.EVT_MENU, self.AddNewHeCell, add_cell)
        #del_cell = self.cellmenu.Append(wx.ID_ANY, '&Delete cell from list', 'Delete cell from list')
        #self.Bind(wx.EVT_MENU, self.DelCell, del_cell)
        
        menubar.Append(filemenu, '&File')
        menubar.Append(self.cellmenu, '&He3 Cells')
        self.SetMenuBar(menubar)
        self.Centre()
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.splitter1 = wx.SplitterWindow(self, -1, style=wx.SP_3D)
        self.splitter2 = wx.SplitterWindow(self.splitter1, -1, style=wx.SP_3D)
        self.splitter3 = wx.SplitterWindow(self.splitter1, -1, style=wx.SP_3D)
        #self.lc1 = wx.ListCtrl(splitter2, -1, style=wx.LC_LIST)
        row = 0
        extent = [ 0, 608, file_objects[row].sample.angle_x.min(), file_objects[row].sample.angle_x.max() ]
        self.fig = Figure(dpi=80, figsize=(5,5))
        self.fig.add_subplot(111)
        self.canvas = Canvas(self.splitter2, -1, self.fig)
        self.canvas.SetMinSize((400,400))
        self.upper_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.lower_sizer = wx.BoxSizer(wx.HORIZONTAL)
        #self.lower_sizer.Add(self.canvas, 1, wx.SHAPED, 0)
        self.lower_sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        self.dataplot = self.fig.axes[0].imshow(log(file_objects[row].detector.counts+1), extent=extent, cmap = cm.jet, aspect='auto', origin='lower', interpolation='nearest')
        panel1 = wx.Panel(self.splitter3, -1, style=wx.SUNKEN_BORDER)
        panel1.SetBackgroundColour("GREY")

        #-- Add the grid to the panel and make it fit
        sizer_panel = wx.BoxSizer(wx.HORIZONTAL)
        sizer_create = wx.FlexGridSizer(rows=6, cols=2, hgap=1, vgap=3)

        self.twoth_zero = wx.TextCtrl(panel1, 1, size=(100, -1))
        twoth_zero_label = wx.StaticText(panel1, -1, 'Main beam\ncenter pixel:')
        sizer_create.Add(twoth_zero_label, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_create.Add(self.twoth_zero)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        self.twoth_zero.SetValue('309.0')

        self.a4_offset = wx.TextCtrl(panel1, 1, size=(100, -1))
        a4_offset_label = wx.StaticText(panel1, -1, 'when A4 = ')
        sizer_create.Add(a4_offset_label, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_create.Add(self.a4_offset)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        self.a4_offset.SetValue('0.0')

        self.bkg_offset = wx.TextCtrl(panel1, 1, size=(100, -1))
        bkg_offset_label = wx.StaticText(panel1, -1, 'Constant bkg.\nsub: ')
        sizer_create.Add(bkg_offset_label, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_create.Add(self.bkg_offset)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        self.bkg_offset.SetValue('0.0')

        self.obj_name = wx.TextCtrl(panel1, 1, size=(100, -1))
        #obj_name_label = wx.StaticText(panel1, -1, 'Object name')
        #sizer_create.Add(obj_name_label, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 0)
        #sizer_create.Add(self.obj_name)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        self.obj_name.SetValue('name')

        self.do_pol_corr = wx.CheckBox(panel1, 1)
        do_pol_corr_label = wx.StaticText(panel1, -1, 'Do Polarization\nCorrection?')
        sizer_create.Add(do_pol_corr_label, 0, wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_create.Add(self.do_pol_corr)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        self.do_pol_corr.SetValue(False)
        #self.Bind(wx.EVT_CHECKBOX, self.toggle_pol_assumptions_active, self.do_pol_corr)
        
        pol_assumption_choices = [
        "no assumptions\n(use all I++, I+-, I-+, I--)",
        "R+- assumed equal to R-+\n(use I++, I-+ and I--)",
        "R-+ assumed equal to R+-\n(use I++, I+- and I--)",
        "R-+ and R+- equal zero\n(use I++, I--)"]
        #self.pol_assumptions = wx.Choice(panel1, 1, choices=pol_assumption_choices)
        #sizer_create.Add(self.pol_assumptions)
        #sizer_create.Add(wx.StaticText(panel1, -1, ''))
        #self.pol_.SetSelection(0)
        #self.pol_assumptions.Disable()

        btn_create = wx.Button(panel1, -1, '&Create combined\nth-2th object')
        self.Bind(wx.EVT_BUTTON, self.CreateDataObject, btn_create)
        sizer_create.Add(btn_create, 0)
        #btn_do_polcorr = wx.Button(panel1, -1, 'Run &Polarization\nCorrection')
        #self.Bind(wx.EVT_BUTTON, self.doPolCorr, btn_do_polcorr
        #sizer_creat.Add(btn_do_polcorr, 0)
        
        sizer_create.Add(self.obj_name)
        
        sizer_panel.Add(sizer_create)
        #sizer_panel.Add(btn_create, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        
        self.obj_tree = wx.TreeCtrl(panel1, 1, wx.DefaultPosition, (-1,-1), wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)
        self.obj_tree.AddRoot('Object tree root')
        self.obj_tree.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnTreeSelChanged, id=1)
        self.obj_tree.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnTreeRightClick)
        sizer_panel.Add(self.obj_tree, -1,  wx.TOP | wx.LEFT | wx.EXPAND)
        panel1.SetSizer(sizer_panel)
        self.panel1 = panel1
        
        #self.IPython_panel = ipython_view.IPShellWidget(self.splitter2, "")
        #self.IPython_panel.SetMinSize((-1,480))
        #self.lower_sizer.Add(self.IPython_panel, 1, wx.RIGHT|wx.TOP)
        blank_panel = wx.Panel(self.splitter2, -1)
        self.lower_sizer.Add(blank_panel, 1, wx.TOP | wx.RIGHT | wx.EXPAND)
        #self.grid = TestGrid(file_objects, file_list, user_descriptions, elapsed_time, self.A3min, self.A3max, self.A4min, self.A4max, parent=self.splitter3, size=(800,300))
        self.grid = GridXML(self.file_catalog, parent=self.splitter3, size=(800,300))
        self.grid.Bind(wx.grid.EVT_GRID_CELL_CHANGE, self.OnGridCellChange)
        self.upper_sizer.Add(self.grid, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        self.upper_sizer.Add(panel1, 0, wx.EXPAND)
        self.splitter3.SetSizer(self.upper_sizer)
        self.splitter2.SetSizer(self.lower_sizer)
        self.splitter2.Fit()
        self.splitter1.SplitHorizontally(self.splitter3, self.splitter2)
        self.splitter3.SplitVertically(self.grid, panel1)
        self.splitter2.SplitVertically(self.canvas, blank_panel)
        self.splitter2.SetSizeHints(400,400,-1,400)
        self.splitter1.SetSashGravity(1.0)
        self.sizer.Add(self.splitter1, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        #splitter2.SplitVertically(dataplot, dataplot)
        #self.Bind(wx.EVT_CLOSE, self.onExit)
        #wx.EVT_CLOSE(self, self.onExit)
        self.fig.canvas.draw()
        self.SetSizer(self.sizer)

        self.prev_rowcol = [None, None]
        wx.EVT_MOTION(self.grid.GetGridWindow(), self.onMotion)
        self.Bind(wx.EVT_CLOSE, close_all)
        # now load up old objects whose state was saved:
        if load_saved_objects:
            self.LoadSavedDataObjects()
    
    def toggle_pol_assumptions_active(self, evt = None):
        if self.do_pol_corr.IsChecked():
            self.pol_assumptions.Enable()
        else:
            self.pol_assumptions.Disable()
     
    def LoadSavedDataObjects(self, filename = None):
        # function to load up all old binned data objects, th2th and qxqz objects, with links to original files
        # throw an error if the original source files are missing.
        # use 'data_objects.xml' to initialize
        if filename == None:
            filename = self.obj_catalog_path
        if not os.access(filename, os.R_OK):
            return # there's no file there, so we can quit looking
        
        data = open(filename, 'r').read() # close is implied
        
        saved_obj_data = xml.dom.minidom.parseString(data)
        obj_catalog = self.obj_catalog.childNodes[0]
        
        cells = saved_obj_data.getElementsByTagName('He3_analyzer')
        for cell in cells:
            #copied_cell = self.obj_catalog.importNode(cell, False)
            #obj_catalog.appendChild(copied_cell)
            
            params = {}
            for key in cell.attributes.keys():
                key = str(key)
                params[key] = cell.getAttribute(key)
            self.AddNewHeCell(params = params)
        
        objs = saved_obj_data.getElementsByTagName('binned_data_object')
        for obj in objs:
            copied_obj = self.obj_catalog.importNode(obj, False)

            file_nodes = obj.getElementsByTagName('binned_data_file')
            files = []
            files_fullpath = []
            file_objects = []
            rows = []
            for file_node in file_nodes:
                copied_file_node = self.obj_catalog.importNode(file_node, False)

                filename = file_node.getAttribute('filename')
                files.append(filename)
                files_fullpath.append(os.path.join(self.data_directory, filename))
                # which row of grid does this correspond to?
                grid_row = None
                for row in range(self.grid.GetNumberRows()):
                    #print self.grid.GetCellValue(row, 0), filename
                    if self.grid.GetCellValue(row, 0) == filename:
                        grid_row = row
                        break
                if grid_row == None:
                    print "Error! Attempt to load data file (%s) that is not in this directory!" % filename
                    return
                else:
                    rows.append(grid_row)
                    copied_obj.appendChild(copied_file_node)
                    file_objects.append(self.file_objects[grid_row])
            print len(file_objects)
            tth0 = float(obj.getAttribute('twoth_zero_pixel'))
            a40 = float(obj.getAttribute('a4_offset'))
            desc = obj.getAttribute('name')
            do_pol_corr = (obj.getAttribute('do_pol_corr') == 'True')
            pol_assum_str = obj.getAttribute('assumptions')
            if pol_assum_str == '':
                pol_assumptions = 0
            else:
                pol_assumptions = int(pol_assum_str)
            new_obj = rebinned_data(file_objects = file_objects, do_pol_corr = do_pol_corr, assumptions = pol_assumptions, twoth_zero_pixel = tth0, a4_offset = a40, description = desc, supervisor = self)
            new_obj._browser_selected_rows = rows
            obj_catalog.appendChild(copied_obj)
            self.SaveObjCatalog()
        
        
            

    def SaveObjCatalog(self, filename = None):
        if filename == None:
            filename = self.obj_catalog_path
        f_out = open(filename, 'w')
        f_out.write(self.obj_catalog.toprettyxml())
        f_out.close()
        
    def onExit(self, event=None):
        self.Destroy()
        
    def AddNewHeCell(self, event=None, params = {}):
        default_params = He3Analyzer.default_params.copy()
        default_params.setdefault('id_num', '1')
        default_params.setdefault('name', '')
        # I'm letting the He3Analyzer class define most of the params, but make sure it has
        # at least an id_num and name for our purposes!
             
        # check to see all needed values are present     
        if not(set(default_params.keys()).issubset(set(params.keys()))):
            default_params.update(params) #get any user-supplied values into default_params
            params = default_params # then use our updated dict as the params
            dlg = get_cell_params_gui2(None, -1, 'Enter parameters for a new cell', params)
            if not dlg.ShowModal() == wx.ID_OK:
                dlg.Destroy()
                return #cancel button pushed
            dlg.Destroy()
            
        # validation of unique id:
        while self.He3_cells.has_key(params['id_num']):
            dlg = get_cell_params_gui2(None, -1, 'Non-unique ID: re-enter parameters for a new cell', params)
            if not dlg.ShowModal() == wx.ID_OK:
                dlg.Destroy()
                return #cancel button pushed
            dlg.Destroy()
        
        if not params == {}:
            cell_id = params['id_num']
            new_cell = He3Analyzer(params)
            self.He3_cells[cell_id] = new_cell
            keys = ['name']
            menustr = '('+cell_id+') '
            for key in keys:
                menustr += key + ': ' + params[key] + ', '
            new_cell_menu = wx.Menu()
            new_cell_menuitem = self.cellmenu.AppendMenu(wx.ID_ANY, menustr, new_cell_menu)
            new_cell_del = new_cell_menu.Append(wx.ID_ANY, '&Delete this cell', 'Delete this cell')
            new_cell_mod = new_cell_menu.Append(wx.ID_ANY, '&Modify this cell', 'Modify cell')
            self.Bind(wx.EVT_MENU, lambda x: self.DelCell(new_cell_menuitem, cell_id), new_cell_del)
            self.Bind(wx.EVT_MENU, lambda x: self.ModifyHeCell(new_cell_menuitem, params), new_cell_mod)
            doc = self.obj_catalog
            new_child = doc.createElement('He3_analyzer')
            for key in params.keys():
                new_child.setAttribute(key, params[key])
            doc.childNodes[0].appendChild(new_child)
            self.SaveObjCatalog()
                        
        
    def ModifyHeCell(self, cell_menuitem, params, event=None):
        #params = self.He3_cells[cell_id]
        dlg = get_cell_params_gui2(None, -1, 'Enter parameters for a new cell', params)
        if not dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy()
            return # do this if cancel button pushed
        dlg.Destroy()
        if not params == {}:
            cell_id = params['id_num']
            new_cell = He3Analyzer(params)
            self.He3_cells[cell_id] = new_cell
            keys = ['name']
            menustr = '('+cell_id+') '
            for key in keys:
                menustr += key + ': ' + params[key] + ', '
            cell_menuitem.SetText(menustr)
            
            doc = self.obj_catalog
            for item in doc.getElementsByTagName('He3_analyzer'):
                if item.getAttribute('id_num') == cell_id:
                    for key in params.keys():
                        item.setAttribute(key, params[key])
            self.SaveObjCatalog()
            
            
        
    def DelCell(self, cell_menuitem, cell_id, event=None):
        self.cellmenu.DeleteItem(cell_menuitem)
        self.He3_cells.pop(cell_id)

    def OnGridCellChange(self, evt):
        print "cell changed: I will save now."
        #figure out which cell changed
        row = evt.GetRow()
        col = evt.GetCol()
        new_data = self.grid.GetCellValue(row,col)
        print "new data: %s" % new_data
        catalog = self.file_catalog.getElementsByTagName('file_catalog')[0]
        for item in catalog.childNodes:
            if item.getAttribute('index') == str(row): #find the appropriate element
                if col == 8:
                    item.setAttribute('user_description', new_data)
                elif col == 1:
                    item.setAttribute('polarization', new_data)
                elif col == 2:
                    item.setAttribute('cell_id', new_data)
                break # we're done with the for loop
        self.SaveFileCatalog()
        
    def SaveFileCatalog(self, filename = None):
        if filename == None:
            filename = self.file_catalog_path
        f_out = open(filename, 'w')
        f_out.write(self.file_catalog.toprettyxml())
        f_out.close()
        
    def loadCatalog(self, filename, rootname = 'file_catalog'):
        f_in = open(filename, 'r')
        data = f_in.read()
        f_in.close()
        old_doc = xml.dom.minidom.parseString(data)
        catalog = old_doc.getElementsByTagName(rootname)[0]
        return catalog
        
    def SaveCatalogSpreadsheet(self, event = None, filename = None, delimiter = ','):
        if filename == None:
            filename = wx.SaveFileSelector('comma-separated value spreadsheet', '.csv', 'file_catalog.csv')
        f_out = open(filename, 'w')
        elements = ['filename', 'start_time', 'polarization', 'cell_id', 'A3min', 'A3max', 'A4min', 'A4max', 'ICP_description', 'user_description']
        element_types = ['string', 'string', 'string', 'int', 'float', 'float', 'float', 'float', 'string', 'string']
        #f_out.write('#') # put in comment prefix
        for element in elements:
            f_out.write(element+delimiter)
        f_out.write('\n')
        for item in self.file_catalog.childNodes[0].childNodes:
            outstr = ''
            for element, etype in zip(elements, element_types):
                if etype == 'string':
                    outstr += '"'+item.getAttribute(element)+'"' + delimiter
                else:
                    outstr += item.getAttribute(element) + delimiter
            outstr += '\n'
            f_out.write(outstr)
        f_out.close()
    
    def GetTimestamps(self, event = None):
        data_dir = self.data_directory
        doc = self.file_catalog
        for item in doc.getElementsByTagName('binned_data_file'):
            filename = item.getAttribute('filename')
            base_file = os.path.join(data_dir, filename[-12:]) # compare to 8.3-format filename without prefixes (I, PI)
            if os.access(base_file, os.R_OK): # if the original file is in this directory, grab the timestamp
                    (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(base_file)
                    start_timestamp = time.mktime(time.strptime(item.getAttribute('start_time')))
                    end_timestamp = mtime # from the os.stat operation above - UNIX/Linux only
                    elapsed = str(end_timestamp - start_timestamp)
                    item.setAttribute('elapsed_time', elapsed)
        self.grid.updateFromXML(self.file_catalog)
        self.SaveFileCatalog()

    def ProcessDirectory(self):
        file_objects = []
        user_descriptions = []
        elapsed_time = []
        A3min = []
        A3max = []
        A4min = []
        A4max = []
        data_dir = self.data_directory
        doc = self.file_catalog
        if os.access(self.file_catalog_path, os.R_OK):
            old_catalog = self.loadCatalog(self.file_catalog_path)
            print "catalog loaded"
        else:
            old_catalog = None
        catalog = doc.createElement('file_catalog')
        catalog.setAttribute('directory', data_dir)
        doc.appendChild(catalog)
        print data_dir
        file_list = glob1(data_dir, 'I*.c?1')
        if len(file_list) < 1: 
            print "no files to display... exiting"
            self.onExit()
        sorting_list = []
        for filename in file_list:
            loaded_file = red.load(os.path.join(data_dir, filename))
            sorting_list.append([filename, loaded_file])
        sorting_list.sort(key = lambda entry: entry[0])
        sorted_list = sorted(sorting_list, key = lambda entry: entry[1].date)
        row = 0
        new_file_list = []
        for item in sorted_list:
            print item[0]
            cat_entry = doc.createElement('binned_data_file')
            cat_entry.setAttribute('filename', item[0])
            cat_entry.setAttribute('index', str(row))
            loaded_file = item[1]
            #loaded_file = red.load(os.path.join(data_dir, item))
            new_file_list.append(item[0])
            file_objects.append(loaded_file)
            cat_entry.setAttribute('ICP_description', loaded_file.description)
            # now grab stuff from the xml file with user-added descriptions in it
            user_description = ''
            elapsed = ''
            
            #def calculate_elapsed(filename):
            #    base_file = os.path.join(data_dir, filename[-12:]) # compare to 8.3-format filename without prefixes (I, PI)
            #    if os.access(base_file, os.R_OK): # if the original file is in this directory, grab the timestamp
            #        (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(base_file)
            #        start_timestamp = time.mktime(time.strptime(loaded_file.date.ctime()))
            #        end_timestamp = mtime # from the os.stat operation above - UNIX/Linux only
            #        elapsed = str(end_timestamp - start_timestamp)
            #        return elapsed
            #   else:
            #        return ''
                   
            if old_catalog:
                for old_item in old_catalog.getElementsByTagName('binned_data_file'):
                    if old_item.getAttribute('filename') == item[0]:
                        if old_item.hasAttribute('user_description'):
                            user_desc = old_item.getAttribute('user_description')
                            cat_entry.setAttribute('user_description', user_desc)
                        if old_item.hasAttribute('polarization'):
                            pol = old_item.getAttribute('polarization')
                            cat_entry.setAttribute('polarization', pol)
                            loaded_file.polarization = pol
                        if old_item.hasAttribute('cell_id'):
                            cell_id = old_item.getAttribute('cell_id')
                            cat_entry.setAttribute('cell_id', cell_id)
                            loaded_file.cell_id = cell_id
                        if old_item.hasAttribute('elapsed_time'): # don't overwrite once the elapsed time has been established!
                            elapsed = old_item.getAttribute('elapsed_time')
                            loaded_file.elapsed_time = elapsed
                        #else:
                        #    elapsed = calculate_elapsed(item[0])
                        break
            #else:
            #    elapsed = calculate_elapsed(item[0])
            if not elapsed == '':
                cat_entry.setAttribute('elapsed_time', elapsed)
                    
            user_descriptions.append(user_description)
            elapsed_time.append(elapsed)
            A3 = loaded_file.sample.angle_x
            cat_entry.setAttribute('A3min', str(A3.min()))
            cat_entry.setAttribute('A3max', str(A3.max()))  
            A3min.append(A3.min())
            A3max.append(A3.max())
            A4 = loaded_file.detector.angle_x
            cat_entry.setAttribute('A4min', str(A4.min()))
            cat_entry.setAttribute('A4max', str(A4.max()))  
            A4min.append(A4.min())
            A4max.append(A4.max())
            cat_entry.setAttribute('start_time', loaded_file.date.ctime())
            catalog.appendChild(cat_entry)
            row += 1 
        self.A3min = min(A3min)
        self.A3max = max(A3max)
        if self.A3max == self.A3min: # need some range to properly display
            self.A3max += 1.0 
            self.A3min -= 1.0 
        self.A4min = min(A4min)
        self.A4max = max(A4max)
        if self.A4max == self.A4min: # need some range to properly display        
            self.A4max += 1.0 
            self.A4min -= 1.0
        
        print len(file_objects), len(new_file_list), len(user_descriptions)
        return file_objects, new_file_list, user_descriptions, elapsed_time

    def onMotion(self, evt):
        x, y = self.grid.CalcUnscrolledPosition(evt.GetPosition())
        row = self.grid.YToRow(y)
        col = self.grid.XToCol(x)
#        if (row,col) != self.prev_rowcol and row >= 0 and col >= 0 and evt.LeftIsDown():
        if (row,col) != self.prev_rowcol and row >= 0 and col >= 0:
            self.prev_rowcol[:] = [row,col]
            file_objects = self.file_objects
            extent = [ 0, 608, file_objects[row].sample.angle_x.min(), file_objects[row].sample.angle_x.max() ]
            if self.dataplot: 
              self.dataplot.set_array(log(file_objects[row].detector.counts + 1))
              self.dataplot.set_extent(extent)
              self.fig.axes[0].set_title(file_objects[row].name)
              self.fig.canvas.draw()
        evt.Skip()        

    def OnTreeSelChanged(self, evt):
        item =  evt.GetItem()
        obj = self.obj_tree.GetPyData(item)
        if isinstance(obj, rebinned_data):
            rows = obj._browser_selected_rows
            self.grid.ClearSelection()
            for row in rows:
                self.grid.SelectRow(row, True)
    
    def OnTreeRightClick(self, evt):
        item = evt.GetItem()
        self.item = item
        obj = self.obj_tree.GetPyData(item)
        if isinstance(obj, rebinned_data):
            name = obj.description
            popup = wx.Menu()
            new_ID = wx.NewId()
            popup.Append(new_ID, '&Subtract Background', 'Subtract Background')
            wx.EVT_MENU(self, new_ID, obj.subtractBackground)
            new_ID = wx.NewId()
            popup.Append(new_ID,'&Convert to Qx-Qz', 'Convert to Qx-Qz')
            wx.EVT_MENU(self, new_ID, obj.runConversion)
            if obj.do_pol_corr:
                new_ID = wx.NewId()
                popup.Append(new_ID, 'Run &Polarization Correction', 'Run Polarization Correction')
                wx.EVT_MENU(self, new_ID, obj.runPolCorr)
            new_ID = wx.NewId()
            popup.Append(new_ID, '&Delete', 'Delete')
            wx.EVT_MENU(self, new_ID, lambda evt: self.DelRebinnedData(name))
            (windowx, windowy) = self.panel1.GetPosition() + self.obj_tree.GetPosition()
            (x, y) = evt.GetPoint() + (windowx, windowy)
            self.PopupMenu(popup, (x, y))
        elif isinstance(obj, plottable_2d_data):
            popup = wx.Menu()
            item1 = popup.Append(wx.ID_ANY,'&Delete', 'Delete')
            wx.EVT_MENU(self, item1.GetId(), lambda evt: self.DelPlottable2dData(obj.number))
            item2 = popup.Append(wx.ID_ANY, '&Replot', 'Replot')
            wx.EVT_MENU(self, item2.GetId(), lambda evt: obj.wxplot() )
            (windowx, windowy) = self.panel1.GetPosition() + self.obj_tree.GetPosition()
            (x, y) = evt.GetPoint() + (windowx, windowy)
            self.PopupMenu(popup, (x, y))
    
    def CreateDataObject(self, evt):
        file_objects = self.file_objects
        rows = self.grid.GetSelectedDataRows()
        files = []
        output_objs = []
        for row in rows:
            print(self.grid.GetCellValue(row, 0))
            files.append(self.grid.GetCellValue(row, 0))
            file_obj = file_objects[row]
            file_obj.elapsed_time = self.grid.GetCellValue(row, 4)
            file_obj.cell_id = self.grid.GetCellValue(row, 2)
            file_obj.polarization = self.grid.GetCellValue(row, 1)
            output_objs.append(file_obj)
                        
        name = self.obj_name.GetValue()
        while self.data_objects.has_key(name):
            dlg = wx.TextEntryDialog(None, "That name already exists.  Enter another:", '')
            if dlg.ShowModal() == wx.ID_OK:
                name = dlg.GetValue()
            dlg.Destroy()
        do_pol_corr = self.do_pol_corr.GetValue()
        #pol_assumptions = self.pol_assumptions.GetSelection()

        #if pol_assumptions == 'no assumptions\n(use all 4 I++, I+-, I-+, I--)':
        #    pol_corr_mode = '++,-+,+-,--'
        #elif pol_assumptions == 'R+- assumed equal to R-+\n(use I++, I-+ and I--)':
        #    pol_corr_mode = '++,-+,--'
        #elif pol_assumptions == 'R-+ assumed equal to R+-\n(use I++, I+- and I--)':
        #    pol_corr_mode = '++,+-,--'
        #elif pol_assumptions == 'R-+ and R+- equal zero\n(use I++, I--)':
        #    pol_corr_mode = '++,--'
        
        # explicitly setting assumptions to None here: change back to self.pol_assumptions.GetSelection() to choose manually.
        new_obj = rebinned_data(file_objects = output_objs, do_pol_corr = do_pol_corr, assumptions = None, twoth_zero_pixel = float(self.twoth_zero.GetValue()), a4_offset = float(self.a4_offset.GetValue()), description = name, plot_data = False, supervisor = self)
        new_obj._browser_selected_rows = rows 
        #self.dataobjects[name] = new_obj
        # add to xml tree of data objects
        doc = self.obj_catalog
        new_child = doc.createElement('binned_data_object')
        new_child.setAttribute('a4_offset', self.a4_offset.GetValue())
        new_child.setAttribute('twoth_zero_pixel', self.twoth_zero.GetValue())
        new_child.setAttribute('name', name)
        new_child.setAttribute('do_pol_corr', str(self.do_pol_corr.GetValue()))
        #new_child.setAttribute('assumptions', str(self.pol_assumptions.GetSelection()))
        doc.childNodes[0].appendChild(new_child)
        for filename in files:
            for item in self.file_catalog.getElementsByTagName('binned_data_file'):
                if item.getAttribute('filename') == filename:
                    #new_child_file = doc.importNode(item, False)
                    new_child_file = doc.createElement('binned_data_file')
                    new_child_file.setAttribute('filename', filename)
                    new_child.appendChild(new_child_file)
        self.SaveObjCatalog()
        #print files
      
    def __getitem__(self, name):
        return self.data_objects[name]
        
    def AddRebinnedData(self, new_object, name = ''):
        """when new object is created, it adds itself to the supervisor list... and this frame acts as a supervisor"""
        self.data_objects[name] = new_object
        new_object.number = self.rb_count
        tree_root = self.obj_tree.GetRootItem()
        branch = self.obj_tree.AppendItem(tree_root, name)
        self.rebinned_data_treeitem[name] = branch       
        self.obj_tree.SetPyData(branch, new_object)
        self.rb_count += 1

        
    def AddPlottable2dData(self, new_object, base_data_obj = None, name=''):
        self.plottable_2d_data_objects.append(new_object)
        new_object.number = self.plottable_count
        tree_root = self.obj_tree.GetRootItem()
        if base_data_obj: tree_root = self.rebinned_data_treeitem[base_data_obj.description]
        subbranch = self.obj_tree.AppendItem(tree_root, name)
        self.plottable_2d_data_treeitem[self.plottable_count] = subbranch
        self.obj_tree.SetPyData(subbranch, new_object)
        self.plottable_count += 1

    def DelPlottable2dData(self, number):
        treeItem = self.plottable_2d_data_treeitem.pop(number)
        obj = self.plottable_2d_data_objects.pop(number)
        self.obj_tree.Delete(treeItem)
        obj.__del__()
        print "deleted object " + str(number)
    
    def DelRebinnedData(self, name):
        obj = self.data_objects[name]
        treeItem = self.rebinned_data_treeitem.pop(name)
        self.obj_tree.Delete(treeItem)
        del(obj)
        # now get it out of the xml tree:
        for child_node in self.obj_catalog.childNodes[0].childNodes: # first child node is docroot
            if child_node.getAttribute('name') == name:
                self.obj_catalog.childNodes[0].removeChild(child_node)
                break
        self.SaveObjCatalog()
        
        
    def onSelectRange(self, evt):
      rows = self.grid.SelectedRows
      if len(rows) == 1:
        row = rows[0]
        file_objects = self.file_objects
        extent = [file_objects[row].detector.angle_x.max() + 7.4, file_objects[row].detector.angle_x.min(), file_objects[row].sample.angle_x.min(), file_objects[row].sample.angle_x.max() ]
        if self.dataplot: 
          self.dataplot.fig.im.set_array(log(file_objects[row].detector.counts + 1))
          self.dataplot.fig.im.set_extent(extent)
          self.dataplot.canvas.draw()
        

class GridXML(wx.grid.Grid):
    """ grid class that uses XML as input/output """
    def __init__(self, file_catalog_xml, parent = None, **kwargs):
        wx.grid.Grid.__init__(self, parent, **kwargs)
        self.updateFromXML(file_catalog_xml)
        
    def updateFromXML(self, file_catalog_xml):
        doc = file_catalog_xml
        elements = doc.getElementsByTagName('binned_data_file')
        self.num_rows = len(elements)
        self.CreateGrid(self.num_rows, 9)
        self.SetColLabelValue(0, 'Filename')
        self.SetColLabelValue(1, 'Pol.')
        self.SetColLabelValue(2, 'Cell')
        self.SetColLabelValue(3, 'Start Time')
        self.SetColLabelValue(4, 'Elaps.(s)')
        self.SetColLabelValue(5, 'A3')
        self.SetColLabelValue(6, 'A4')
        self.SetColLabelValue(7, 'ICP_Descriptor')
        self.SetColLabelValue(8, 'Description')
        self.SetRowLabelSize(40)
        self.tooltips = []
        # get max and mins:
        A3mins = []
        A3maxs = []
        A4mins = []
        A4maxs = []
        for item in elements:
            A3mins.append(float(item.getAttribute('A3min')))
            A3maxs.append(float(item.getAttribute('A3max')))
            A4mins.append(float(item.getAttribute('A4min')))
            A4maxs.append(float(item.getAttribute('A4max')))
        A3mins = numpy.array(A3mins)
        A3maxs = numpy.array(A3maxs)
        A4mins = numpy.array(A4mins)
        A4maxs = numpy.array(A4maxs)
        A3min = A3mins.min()
        A3max = A3maxs.max()
        A4min = A4mins.min()
        A4max = A4maxs.max()
        if A4max == A4min:
            A4max += 0.0001
            
        for item in elements:
            row = int(item.getAttribute('index'))
            self.SetCellValue(row, 0, item.getAttribute('filename'))
            self.SetCellValue(row, 1, item.getAttribute('polarization')) # will default to '' if not filled
            self.SetCellValue(row, 2, item.getAttribute('cell_id'))
            self.SetCellValue(row, 3, item.getAttribute('start_time'))
            self.SetCellValue(row, 4, item.getAttribute('elapsed_time'))
            self.SetCellValue(row, 7, item.getAttribute('ICP_description'))
            self.SetCellValue(row, 8, item.getAttribute('user_description'))

            self.SetCellEditor(row, 1, wx.grid.GridCellChoiceEditor(['', '++', '-+', '--', '+-']))
            self.SetCellEditor(row, 2, wx.grid.GridCellNumberEditor())

            imgA3 = GetBarGraphBitmap(A3min, A3max, float(item.getAttribute('A3min')), float(item.getAttribute('A3max')), color_on = wx.RED)
            imageRendererA3 = MyImageRenderer(imgA3)
            self.SetCellRenderer(row, 5, imageRendererA3)
            self.SetColSize(5,imgA3.GetWidth()+2)
            #self.grid.SetRowSize(row,img.GetHeight()+2)
            #self.grid.SetCellValue(row, 1, str(file_objects[row].sample.angle_x.min()) + ' to ' + str(file_objects[row].sample.angle_x.max()) )
            A3tooltip = item.getAttribute('A3min') + ' to ' + item.getAttribute('A3max')

            imgA4 = GetBarGraphBitmap(A4min, A4max, float(item.getAttribute('A4min')), float(item.getAttribute('A4max')), color_on = wx.BLUE)
            imageRendererA4 = MyImageRenderer(imgA4)
            self.SetCellRenderer(row, 6, imageRendererA4)
            self.SetColSize(6,imgA4.GetWidth()+2)
            A4tooltip = item.getAttribute('A4min') + ' to ' + item.getAttribute('A4max')
            self.tooltips.append([item.getAttribute('filename'), '', '', item.getAttribute('start_time'),item.getAttribute('elapsed_time'), A3tooltip, A4tooltip, item.getAttribute('ICP_description'), item.getAttribute('user_description')])
            self.SetReadOnly(row, 0, True)
            self.SetReadOnly(row, 3, True)
            self.SetReadOnly(row, 4, True)
            self.SetReadOnly(row, 5, True)
            self.SetReadOnly(row, 6, True)
            self.SetReadOnly(row, 7, True)
        self.AutoSizeColumn(0)
        self.SetColSize(1, 50)
        self.SetColSize(2, 50)
        self.AutoSizeColumn(3)
        self.AutoSizeColumn(4)
        self.AutoSizeColumn(7)
        self.SetColSize(8, 250)

    def GetSelectedDataRows(self):
        """ need to add up all the different ways the selection is stored: in rows, columns,
            individual cells and blocks of cells """
        selection = set([])
        if self.SelectedCols != []:
          return range(self.num_rows)
        else:
          selection = selection.union( self.GetSelectedRows())
          for cell in self.SelectedCells:
            selection.add(cell[0])
          if self.GetSelectionBlockTopLeft() != []:
            selection = selection.union( range(self.GetSelectionBlockTopLeft()[0][0], self.GetSelectionBlockBottomRight()[0][0] + 1) )
          return list(selection)
         
class get_cell_params_gui2(wx.Dialog):
    def __init__(self, parent, id, title, params):
        wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,400), style = wx.DEFAULT_DIALOG_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)
        
        #-- I think you need a panel if you're going to have more than one control in your frame.
        #panel = wx.Panel(self, -1)
        self.params = params

        #-- Create the processing button, add it to the panel and wire it up to a function in the class
        btn_SaveExit = wx.Button(self, wx.ID_OK, "&Done")
        self.Bind(wx.EVT_BUTTON, self.saveExit, btn_SaveExit)
        
        #-- Create the close button, add it to the panel and wire it up to a function in the class
        btn_Close = wx.Button(self, -1, "&Cancel")
        self.Bind(wx.EVT_BUTTON, self.onExit, btn_Close)
        
        #-- Now we have to create a grid to layout our controls
        sizer_main = wx.BoxSizer(wx.VERTICAL)
        #rows=4,cols=1,hgap=1,vgap=5)
        
        sizer_buttons = wx.FlexGridSizer(rows=1,cols=2, hgap=5,vgap=5)
        sizer_buttons.Add(btn_SaveExit)
        sizer_buttons.Add(btn_Close)
        
        sizer_params = wx.FlexGridSizer(rows=20, cols=2, hgap=5, vgap=10)
        self.values = {}
        self.text_labels = {}
        
        self.keys = [
        'P0',
        'T0',
        'Gamma',
        'sigmanl',
        'Psm',
        'Pf',
        'name',
        'id_num',
        ]
        
        self.labels = [
        'Initial Polarization (P0): ',
        'Initial Transmission (T0): ',
        'Decay time const. (Gamma, hours): ',
        'Sigma n l : ', 
        'Supermirror Polarization (Psm): ',
        'Flipper polarization (Pf): ',
        'Cell name (for reference only): ',
        'Cell id# (integer):',
        ]
        for key,label in zip(self.keys,self.labels):
            value = wx.TextCtrl(self, 1, size=(100, -1))
            text_label = wx.StaticText(self, -1, label)
            self.values[key] = value 
            self.text_labels[key] = text_label
            sizer_params.Add(text_label)
            sizer_params.Add(value)
            value.SetValue(str(params[key]))
        
        #self.filenames = params['inFileNames']
        self.static_line_1 = wx.StaticLine(self, -1)
        self.time_label = wx.StaticText(self, -1, "Select t0 (cell initial characterization)")
        self.calendar_ctrl_1 = wx.calendar.CalendarCtrl(self, -1)
        self.spin_ctrl_3 = wx.SpinCtrl(self, -1, "", min=0, max=23)
        self.label_1 = wx.StaticText(self, -1, "Hour")
        self.spin_ctrl_4 = wx.SpinCtrl(self, -1, "", min=0, max=60)
        self.label_2 = wx.StaticText(self, -1, "Minute")
        
        try: 
            t0 = datetime.datetime.strptime(params['t0_str'], '%a %b %d %H:%M:%S %Y')
        except:
            t0 = datetime.datetime.now()
            
        self.spin_ctrl_3.SetValue(t0.hour)
        self.spin_ctrl_4.SetValue(t0.minute)
        wxDate = wx.DateTime()
        wxDate.Set(year = t0.year, month = t0.month -1, day = t0.day)
        self.calendar_ctrl_1.SetDate(wxDate)
        
        
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_3 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5 = wx.BoxSizer(wx.VERTICAL)
        sizer_4 = wx.BoxSizer(wx.VERTICAL)
        sizer_4.Add(self.spin_ctrl_3)
        sizer_4.Add(self.label_1)
        sizer_3.Add(sizer_4)
        sizer_5.Add(self.spin_ctrl_4)
        sizer_5.Add(self.label_2)
        sizer_3.Add(sizer_5)
        sizer_2.Add(self.static_line_1, flag = wx.EXPAND)
        sizer_2.Add(self.time_label)
        sizer_2.Add(sizer_3)
        sizer_2.Add(self.calendar_ctrl_1)
     
        sizer_main.Add(sizer_params)
        sizer_main.Add(sizer_2)
        sizer_main.Add(sizer_buttons)

        self.SetSizer(sizer_main)
        sizer_main.Fit(self)

        #-- Show the window that we've just built    
        
    def saveExit(self, event):
        for key in self.keys:
            self.params[key] = self.values[key].GetValue()
        wxDate = self.calendar_ctrl_1.GetDate() # have to convert to datetime
        t0 = datetime.datetime(wxDate.GetYear(), wxDate.GetMonth() +1, wxDate.GetDay(), self.spin_ctrl_3.GetValue(), self.spin_ctrl_4.GetValue())
        self.params['t0_str'] = t0.ctime()
        self.EndModal(wx.ID_OK)
        #self.Close(True)
      
    def onExit(self, event):
        self.Close(True) 
     
def close_all(evt=None):
    for window in wx.GetTopLevelWindows():
        names = ['DataBrowser', 'offspec_plot_2d_data']
        if window.__class__.__name__ in names:
            window.Destroy()
        

class MyImageRenderer(wx.grid.PyGridCellRenderer):
     def __init__(self, img):
         wx.grid.PyGridCellRenderer.__init__(self)
         self.img = img
     def Draw(self, grid, attr, dc, rect, row, col, isSelected):
         image = wx.MemoryDC()
         image.SelectObject(self.img)
         dc.SetBackgroundMode(wx.SOLID)
         if isSelected:
             dc.SetBrush(wx.Brush(wx.BLUE, wx.SOLID))
             dc.SetPen(wx.Pen(wx.BLUE, 1, wx.SOLID))
         else:
             dc.SetBrush(wx.Brush(wx.WHITE, wx.SOLID))
             dc.SetPen(wx.Pen(wx.WHITE, 1, wx.SOLID))
         dc.DrawRectangleRect(rect)
         width, height = self.img.GetWidth(), self.img.GetHeight()
         if width > rect.width-2:
             width = rect.width-2
         if height > rect.height-2:
             height = rect.height-2
         dc.Blit(rect.x+1, rect.y+1, width, height, image, 0, 0, wx.COPY, True)

def InstallGridHint(grid, tooltips_array):
    prev_rowcol = [None,None]
    def OnMouseMotion(evt):
        # evt.GetRow() and evt.GetCol() would be nice to have here,
        # but as this is a mouse event, not a grid event, they are not
        # available and we need to compute them by hand.
        x, y = grid.CalcUnscrolledPosition(evt.GetPosition())
        row = grid.YToRow(y)
        col = grid.XToCol(x)

        if (row,col) != prev_rowcol and row >= 0 and col >= 0:
            prev_rowcol[:] = [row,col]
            #hinttext = rowcolhintcallback(row, col)
            hinttext = tooltips_array[row][col]
            if hinttext is None:
                hinttext = ''
            grid.GetGridWindow().SetToolTipString(hinttext)
        evt.Skip()

    wx.EVT_MOTION(grid.GetGridWindow(), OnMouseMotion)

def InstallAutosaveDescriptions(grid):
    def onDescriptionChange(evt):
        pass

def GetBitmap( width=32, height=32, colour = (0,0,0) ):
  array = numpy.zeros( (height, width, 3),'uint8')
  array[:,:,] = colour
  image = wx.EmptyImage(width,height)
  image.SetData( array.tostring())
  return image.ConvertToBitmap() # wx.BitmapFromImage(image)

def GetBarGraphBitmap( range_min = -20.0, range_max = 20.0, value_min = 0.0, value_max = 8.0, bmp_width = 50, bmp_height = 20, color_on = (255,0,0), color_off = (255,255,255) ):
  array = numpy.zeros( (bmp_height, bmp_width, 3),'uint8')
  array[:,:,] = color_off
  range_min = float(range_min)
  range_max = float(range_max)
  value_min = float(value_min)
  value_max = float(value_max)
  min_index = numpy.int( (value_min - range_min) / (range_max - range_min) * bmp_width ) 
  max_index = numpy.int( (value_max - range_min) / (range_max - range_min) * bmp_width )
 
  if min_index < 0: min_index = 0
  if min_index > bmp_width: min_index = bmp_width
  if max_index < 0: max_index = 0
  if max_index > bmp_width: max_index = bmp_width

  if min_index == max_index:
    if min_index == bmp_width: 
      min_index = bmp_width - 1
    else:
      max_index = min_index + 1

  array[:,min_index:max_index,] = color_on
  image = wx.EmptyImage(bmp_width,bmp_height)
  image.SetData( array.tostring())
  return image.ConvertToBitmap() # wx.BitmapFromImage(image)

from ConfigParser import ConfigParser
Config = ConfigParser()
#

def new_databrowser(data_directory = None):
    Config = ConfigParser()
    if __name__ == '__main__':
        config_dir = os.path.split(sys.argv[0])[0]
        config_path = os.path.join(config_dir, 'config.ini')
    
    #Else the module was imported and it has a __file__ attribute that will be the full path of the module.
    else:
        config_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    
    #config_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    if os.access(config_path, os.R_OK):
        Config.read(config_path)
        try:
            default_working_directory = Config.get('Paths', 'default_working_directory')
        except:
            default_working_directory = os.getcwd()
    else:
        default_working_directory = os.getcwd()
    
    if data_directory is None:
        dlg = wx.DirDialog(None, "Choose a working directory (where the rebinned 2d data files are)", defaultPath = default_working_directory )
        if dlg.ShowModal() == wx.ID_OK:
            data_directory = dlg.GetPath()
        dlg.Destroy()
    try:
        config_out = open(config_path, 'w')
        if not ('Paths' in Config.sections()):
            Config.add_section('Paths')
        Config.set('Paths', 'default_working_directory', data_directory )
        Config.write(config_out)
        config_out.close()
    except:
        pass
    db = DataBrowser(data_directory = data_directory)
    global supervisor
    supervisor = db
    InstallGridHint(db.grid, db.grid.tooltips)
    db.Show()
    return db

def md5_for_file(f, block_size=2**20):
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.digest()

    
if __name__ == '__main__':
    app = wx.PySimpleApp(redirect=True)
    #app = wx.App(redirect=False)
    db = new_databrowser()
    app.SetTopWindow(db)
    db.Show()
    app.MainLoop()

