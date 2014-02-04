#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.6.3 on Wed Feb  9 16:00:13 2011

import wx

# begin wxGlade: extracode
# end wxGlade



class BkgDialog(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: MyDialog.__init__
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        
        if 'pol_states_available' in kwds:
            self.pol_states_available = kwds.pop('pol_states_available')
        else:
            self.pol_states_available = ['++', '+-', '-+', '--']
        
        wx.Dialog.__init__(self, *args, **kwds)
          
        self.label_1 = wx.StaticText(self, -1, "Background: \navg. over range of theta and subtract\n")
        #self.label_2 = wx.StaticText(self, -1, "Pol state to get bkg from:")
        #self.pol_state = wx.Choice(self, -1, choices=self.pol_states_available)
        self.label_3 = wx.StaticText(self, -1, "Min theta:")
        self.text_ctrl_1 = wx.TextCtrl(self, -1, "")
        self.label_4 = wx.StaticText(self, -1, "Max theta:")
        self.text_ctrl_2 = wx.TextCtrl(self, -1, "")
        self.label_5 = wx.StaticText(self, -1, "Constant bkg:")
        self.text_ctrl_3 = wx.TextCtrl(self, -1, "")
        self.button_1 = wx.Button(self, -1, "Do subtraction")
        self.Bind(wx.EVT_BUTTON, self.Finished, self.button_1)
        self.button_2 = wx.Button(self, -1, "&Cancel")
        self.Bind(wx.EVT_BUTTON, self.onExit, self.button_2)

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: MyDialog.__set_properties
        self.SetTitle("Background")
        self.text_ctrl_1.SetValue('0.0')
        self.text_ctrl_2.SetValue('0.0')
        self.text_ctrl_3.SetValue('0.0')
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: MyDialog.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.FlexGridSizer(4, 2, 10, 10)
        #grid_sizer_1.Add(self.label_2, 0, 0, 0)
        #grid_sizer_1.Add(self.pol_state, 0, 0, 0)
        grid_sizer_1.Add(self.label_3, 0, 0, 0)
        grid_sizer_1.Add(self.text_ctrl_1, 0, 0, 0)
        grid_sizer_1.Add(self.label_4, 0, 0, 0)
        grid_sizer_1.Add(self.text_ctrl_2, 0, 0, 0)
        grid_sizer_1.Add(self.label_5, 0, 0, 0)
        grid_sizer_1.Add(self.text_ctrl_3, 0, 0, 0)
        grid_sizer_1.Add(self.button_2, 0, 0, 0)
        grid_sizer_1.Add(self.button_1, 0, 0, 0)
        sizer_1.Add(self.label_1, 0, 0, 0)
        sizer_1.Add(grid_sizer_1, 1, wx.EXPAND, 0)
        #sizer_1.Add(self.button_1, 0, 0, 0)
        #sizer_1.Add(self.button_2, 0, 0, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade

    def Finished(self, evt = None):
        self.min_theta = float(self.text_ctrl_1.GetValue())
        self.max_theta = float(self.text_ctrl_2.GetValue())
        self.const_bkg = float(self.text_ctrl_3.GetValue())
        #self.pol_state_chosen = self.pol_states_available[self.pol_state.GetSelection()]
        self.EndModal(wx.ID_OK)
        
    def onExit(self, event):
        self.EndModal(wx.ID_CANCEL)
        self.Close(True) 
# end of class MyDialog


if __name__ == "__main__":
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    dialog_1 = BkgDialog(None, -1, "Background: \navg. over range of theta and subtract")
    app.SetTopWindow(dialog_1)
    dialog_1.Show()
    app.MainLoop()
