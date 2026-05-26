#!/usr/bin/env python3

"""
GUI for data file and test platform input for the HopPyBar python program.
TODO: Updated from PyQt5 to PyQt6 on Jul 31, 2025; Still considering replacing everything with PySide of something else
@author: morrow
Created 2020/09/24
"""

try: # Allow for both PyQt6 (default conda channel) and PyQt5 (conda-forge or WinPython)
    from PyQt6.QtCore import (Qt, QSettings)
    from PyQt6.QtWidgets import (QApplication, QDialog, QTextEdit, QFileDialog,
                                QLabel, QHBoxLayout, QPushButton, QGroupBox,
                                QCheckBox, QComboBox, QVBoxLayout, QGridLayout, QMessageBox,
                                QRadioButton, QGridLayout)
except:
    from PyQt5.QtCore import (Qt, QSettings)
    from PyQt5.QtWidgets import (QApplication, QDialog, QTextEdit, QFileDialog,
                                QLabel, QHBoxLayout, QPushButton, QGroupBox,
                                QCheckBox, QComboBox, QVBoxLayout, QGridLayout, QMessageBox,
                                QRadioButton, QGridLayout)
import sys, os
from datetime import datetime

now = datetime.now()
date = '{:02d}{:02d}{:02d}'.format(now.year, now.month, now.day)

class ImportData(QDialog):
    def __init__(self, parent=None):
        super(ImportData, self).__init__(parent)
#
        self.settings = QSettings('HopPyBar', 'IO_Select')
        self.setStyleSheet('QWidget {font: "Arial"}')
        self.initUI()
        self.loadSettings()
#
    def initUI(self):
        self.setWindowTitle("HopPyBar Data Input")
        self.createTopLeftGroupBox()
        self.createTopCenterGroupBox()
        self.createDispBox()
        self.createOutOptionsGroupBox()
        self.createBottomGroupBox()
        self.createOKCancel()
#
        bottomLayout = QHBoxLayout()
        bottomLayout.addWidget(self.bottomGroupBox)
#
        mainLayout = QGridLayout()
        #mainlayout.addWidget(name, row, column, rowspan (opt), columnspan (opt))
        mainLayout.addWidget(self.topLeftGroupBox, 0, 0, 2, 1) # spanning 2 rows (one column)
        mainLayout.addWidget(self.topCenterGroupBox, 0, 1)
        mainLayout.addWidget(self.dispBox, 0, 2)
        mainLayout.addWidget(self.outOptionsGroupBox, 1, 2)
        mainLayout.addLayout(bottomLayout, 2, 0, 1, 2)
        mainLayout.addWidget(self.OKCancel, 2, 2)
        #mainLayout.setRowStretch(1, 1)
        #mainLayout.setRowStretch(2, 1)
        mainLayout.setColumnStretch(0, 1)
        mainLayout.setColumnStretch(1, 1)
        self.setLayout(mainLayout)
#
    def loadSettings(self):
        try:
            self.dispersion_select.setCurrentIndex(self.settings.value('default_dispersion_select'))
        except:
            pass
        try:
            self.suppress_output_check.setChecked(self.settings.value('default_suppress_output', type=bool))
            self.suppress_granta_check.setChecked(self.settings.value('default_suppress_granta', type=bool))
            self.eng_out_check.setChecked(self.settings.value('default_eng_out', type=bool))
        except:
            self.suppress_output_check.setChecked(False)
            self.suppress_granta_check.setChecked(False)
            self.eng_out_check.setChecked(False)   
#
    def openFileNameDialog(self):
        self.fileid, _ = QFileDialog.getOpenFileName(self,
                            "Select Input File",
#                            "",
                            self.settings.value('recent_dir'),
                            "All Files (*);;Raw Files (*.raw);;Text Files (*.txt);;CSV Files (*.csv)",
                            options=QFileDialog.Option.DontUseNativeDialog)
        #self.filename.setText(fileid)
        if self.fileid:
            self.filename.setText(self.fileid)
            self.data_dir, filename = os.path.split(self.fileid)  # same thing as dirname, but also gets us filename for free
            dirlist = os.listdir(self.data_dir)
            if not os.path.exists(self.data_dir + os.sep + 'Analysis_{}'.format(str(date))):
                self.overwrite = None
                self.bottomGroupBox.hide()  # hides the overwrite dialog if we don't need it, in case the user canceled the file selection after a 2nd attempt.
            else:
                self.bottomGroupBox.show()
#                self.OverwriteDialog()
        else:
            print('File selection cancelled.')
            self.fileid = None
            self.overwrite = None
    def createTopLeftGroupBox(self):
        self.topLeftGroupBox = QGroupBox("Select Input")
        filetext = QLabel('Input File:')
        SelectFilePushButton = QPushButton("Select")
        SelectFilePushButton.setDefault(False)
        self.filename = QTextEdit("File not selected.\n"
                                "(Note: Select any channel; "
                                "other channels are imported automatically)")
        SelectFilePushButton.clicked.connect(self.openFileNameDialog)
#
        layout = QVBoxLayout()
        layout.addWidget(filetext)
        layout.addWidget(self.filename)
        layout.addWidget(SelectFilePushButton)
        layout.addStretch(1)
        self.topLeftGroupBox.setLayout(layout)
#
    def createTopCenterGroupBox(self):
        self.topCenterGroupBox = QGroupBox("Select Test Platform")
        self.topCenterGroupBox.setCheckable(False)
        self.topCenterGroupBox.setFlat(True)
#
        vbox = QVBoxLayout()
        self.topCenterGroupBox.setLayout(vbox)
#
        radiobutton = QRadioButton("LANL Kolsky Systems")
        radiobutton.platform = "LANL"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)
#
        radiobutton = QRadioButton("REL Kolsky Systems")
        radiobutton.platform = "REL"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)
#
        radiobutton = QRadioButton("Mini-Kolsky")
        radiobutton.setChecked(False)
        radiobutton.platform = "mini"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)   #layout.addWidget(radiobutton, 0, 0)
#
        radiobutton = QRadioButton("LANL FLAG Output")
        radiobutton.platform = "FLAG"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)
 #
        radiobutton = QRadioButton("APS/DCS")
        radiobutton.platform = "APS"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)
#
        radiobutton = QRadioButton("PDV Data")
        radiobutton.platform = "PDV"
        if self.settings.value('default_platform') == radiobutton.platform:
            radiobutton.setChecked(True)
            self.platform = radiobutton.platform
        radiobutton.toggled.connect(self.onClicked)
        vbox.addWidget(radiobutton)
#
    def onClicked(self):
        radioButton = self.sender()
        if radioButton.isChecked():
            self.platform = radioButton.platform
#
    def createDispBox(self):
        dbox = QGridLayout()
        self.dispBox = QGroupBox("Dispersion Correction")
        self.dispBox.setLayout(dbox)
#
        text_ds = QLabel("Method")
        self.dispersion_select = QComboBox(self)
        self.dispersion_select.addItem("Bragov 2022")
        self.dispersion_select.addItem("None")
        self.dispersion_select.setCurrentIndex(0)
        dbox.addWidget(text_ds, 1, 0)
        dbox.addWidget(self.dispersion_select, 1, 1)
        #
#    
    def createOutOptionsGroupBox(self):
        self.outOptionsGroupBox = QGroupBox("Output Options")
#
        #vbox = QVBoxLayout()
        vbox = QGridLayout()
        self.outOptionsGroupBox.setLayout(vbox)
#
        self.suppress_output_check = QCheckBox("Suppress verbose output (many files/plots)")
        self.suppress_output_check.setCheckState(Qt.CheckState.Unchecked)   # set default to unchecked
        self.suppress_output_check.stateChanged.connect(self.onChecked)
        vbox.addWidget(self.suppress_output_check, 2, 0)   #layout.addWidget(radiobutton, 0, 0)
#
        self.suppress_granta_check = QCheckBox("Suppress GRANTA Output")
        self.suppress_granta_check.setCheckState(Qt.CheckState.Unchecked)   # set default to unchecked
        self.suppress_granta_check.stateChanged.connect(self.onChecked2)
        vbox.addWidget(self.suppress_granta_check, 3, 0)   #layout.addWidget(radiobutton, 0, 0)
#
        self.eng_out_check = QCheckBox("Output Engineering Values Instead of True")
        self.eng_out_check.setCheckState(Qt.CheckState.Unchecked)   # set default to unchecked
        self.eng_out_check.stateChanged.connect(self.onChecked3)
        vbox.addWidget(self.eng_out_check, 4, 0)   #layout.addWidget(radiobutton, 0, 0)
#

        if self.suppress_output_check.isChecked() == True:
            self.suppress_output = 2
        else:
            self.suppress_output = 0
        if self.suppress_granta_check.isChecked() == True:
            self.suppress_granta = 2
        else:
            self.suppress_granta = 0
        if self.eng_out_check.isChecked() == True:
            self.eng_out = 2
        else:
            self.eng_out = 0
#
    def onChecked(self):
        checkBox = self.sender()
        self.suppress_output = checkBox.checkState()

    def onChecked2(self):
        checkBox = self.sender()
        self.suppress_granta = checkBox.checkState()

    def onChecked3(self):
        checkBox = self.sender()
        self.eng_out = checkBox.checkState()
#
    def createBottomGroupBox(self):
        self.bottomGroupBox = QGroupBox("")
        self.bottomGroupBox.setFlat(True)
        self.bottomGroupBox.setStyleSheet("background-color: red;")
        bbox = QVBoxLayout()
        self.bottomGroupBox.setLayout(bbox)
        self.overwriteText = QLabel('WARNING: Analysis directory already exists!')
        self.overwriteDataCheck = QCheckBox("&Overwrite Analysis Data? (Verbose mode only)")
        self.overwriteDataCheck.setChecked(False)
        self.overwrite = 'no'      #initialize to no overwrite, just incase it never gets clicked and runs the check.
        self.overwriteDataCheck.toggled.connect(self.overwriteClicked)
        bbox.addWidget(self.overwriteText)
        bbox.addWidget(self.overwriteDataCheck)
        #self.bottomGroupBox.setDisabled(True)   # disabled until we know we need it. Gets enabled by overwriteDialog()
        self.bottomGroupBox.hide()
#
    def overwriteClicked(self):
        self.overwriteDataCheck = self.sender()
        if self.overwriteDataCheck.isChecked():
            self.overwrite = 'yes'
        else:
            self.overwrite = 'no'
    def inputErrorDialog(self):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText("Input error, please select file or test platform.")
        msgBox.setWindowTitle("Warning")
        msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        #
        returnValue = msgBox.exec()
        if returnValue == QMessageBox.Ok:
            return
        elif returnValue == QMessageBox.Cancel:
            #print('Aborting data analysis.')
            #sys.exit()
            return
#
    def OKcontinue(self, event):
        if not (hasattr(self,'platform') and hasattr(self, 'data_dir')):#'data_dir' in locals()):   #'platform' in globals() and
            self.inputErrorDialog()
        else:
            self.settings.setValue('default_platform', self.platform)
            self.settings.setValue('default_dispersion_select', self.dispersion_select.currentIndex())
            self.settings.setValue('default_suppress_output', self.suppress_output_check.isChecked())
            self.settings.setValue('default_suppress_granta', self.suppress_granta_check.isChecked())
            self.settings.setValue('default_eng_out', self.eng_out_check.isChecked())
            self.settings.setValue('recent_dir', self.data_dir)
            print('Input selection complete. Default settings saved to: {}'.format(self.settings.fileName()))
            super(ImportData, self).accept()
            #self.close()

    def cancel(self):
        print('Input selection aborted by user.')
        self.close()
        sys.exit()
#
    def createOKCancel(self):
        self.OKCancel = QGroupBox()
        self.OKCancel.setFlat(True)
        OKPushButton = QPushButton("OK (Continue)")
        OKPushButton.setDefault(True)
        OKPushButton.clicked.connect(self.OKcontinue)
        CancelPushButton = QPushButton("Cancel (Quit)")
        CancelPushButton.setDefault(False)
        CancelPushButton.clicked.connect(self.cancel)
#
        layout = QVBoxLayout()
        layout.addWidget(OKPushButton)
        layout.addWidget(CancelPushButton)
        #layout.addStretch(1)
        self.OKCancel.setLayout(layout)
#
    def get_output(self):
            print(f'verbose: {self.suppress_output_check.isChecked()}')
            print(f'granta: {self.suppress_granta_check.isChecked()}')
            print(f'eng: {self.eng_out_check.isChecked()}')
            return (self.fileid, self.platform, self.overwrite, self.dispersion_select, self.suppress_output_check.isChecked(),
                self.suppress_granta_check.isChecked(), self.eng_out_check.isChecked(), self.data_dir)
#
if __name__ == '__main__':
    import os, sys
    app = QApplication(sys.argv)
    import_ui = ImportData()
    import_ui.show()
    sys.exit(app.exec_())
