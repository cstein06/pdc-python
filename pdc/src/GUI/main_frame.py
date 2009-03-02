import wx
import os

from data_io.eeg_in import read_hor_eeg

class MainWindow(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title)
        self.build_menu()
        
    def build_menu(self):
        file_menu = wx.Menu()
        file_menu.Append(ID_OPEN_32H, "&Open 32-C Horizontal EEG",
                         "Open EEG data with a channel per row")
        file_menu.Append(ID_EXIT, "Exit", "Terminate the program")
        
        wx.EVT_MENU(self, ID_OPEN_32H, self.open_eeg_32h)
        wx.EVT_MENU(self, ID_EXIT, self.exit)
        
    def open_eeg_32h(self):
        "Open 32-C EEG data with a channel per row"
        self.dirname = os.getcwd()
        open_dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if (open_dlg.ShowModal() == wx.ID_OK):
            self.filename = open_dlg.GetFilename()
            self.dirname = open_dlg.GetDirectory()
            self.abs_filename = os.path.join(self.filename, self.dirname)
            wait_dialog = wx.MessageDialog