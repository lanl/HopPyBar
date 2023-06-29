#!/usr/bin/env python3

"""
GUI for Granta-specific parameters for the HopPyBar python program.
This enables us to skip the DA reprocessing step.
@author: morrow
Created 2020/12/15
"""
#TODO: make fancier/functional; add ability to do direct imports on Windows/Linux (no support for GRANTA STK toolkit on macOS) 

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import * #QDialog, QInputDialog, QLineEdit

import sys, os
import pandas as pd

#%%########## PULL IN GRANTA OPTIONS FROM LISTS ##########
granta_dir = os.path.join('.','GRANTA') + os.sep 
granta_user_list = pd.read_csv(granta_dir + 'users.txt', delimiter = '\t', header = 0, names = ['user', 'group', 'z_number'])
granta_data_sensitivity = pd.read_csv(granta_dir + 'granta_dataSensitivity.txt', delimiter = '\t', header = 0, names = ['level'])
granta_test_lab = pd.read_csv(granta_dir + 'granta_testLab.txt', delimiter = '\t', header = 0, names = ['lab'])
granta_test_group = pd.read_csv(granta_dir + 'granta_testGroup.txt', delimiter = '\t', header = 0, names = ['group'])
granta_test_location = pd.read_csv(granta_dir + 'granta_testLocation.txt', delimiter = '\t', header = 0, names = ['location'])
granta_specimen_orientation = pd.read_csv(granta_dir + 'granta_specimenOrientation.txt', delimiter = '\t', header = 0, names = ['orientation'])

class InputGrantaData(QDialog):
    def __init__(self, parent=None):
        super(InputGrantaData, self).__init__(parent)
#
        self.settings = QSettings('HopPyBar', 'Granta_Input')

        self.createTopGroupBox()
        self.createOKCancel()
#
        mainLayout = QGridLayout()
        mainLayout.addWidget(self.topGroupBox, 1, 0)
#        mainLayout.addWidget(self.topRightGroupBox, 1, 1)
#        mainLayout.addLayout(bottomLayout, 2, 0, 1, 1)
        mainLayout.addWidget(self.OKCancel, 2, 0)
        #mainLayout.setRowStretch(1, 1)
        #mainLayout.setRowStretch(2, 1)
        mainLayout.setColumnStretch(0, 1)
        mainLayout.setColumnStretch(1, 1)
        self.setLayout(mainLayout)
#
        self.setWindowTitle("Granta Info Input")
#
    def createTopGroupBox(self):
        self.topGroupBox = QGroupBox("Granta Info")

        text_ds = QLabel('Data Sensitivity:')
        self.dataSensitivity = QComboBox(self)
        for x in range(len(granta_data_sensitivity)):
            self.dataSensitivity.addItem(granta_data_sensitivity.level[x])
        
        self.dataSensitivity.setCurrentIndex(0)  # set default to Non-Sensitive, or whatever is first choice in list

        text_do = QLabel('Data Originator:')
        self.dataOriginator = QComboBox(self)
        for x in range(len(granta_user_list)):
            self.dataOriginator.addItem(granta_user_list.user[x])

        text_tc = QLabel('Testing Contact:')
        self.testContact = QComboBox(self)
        for x in range(len(granta_user_list)):
            self.testContact.addItem(granta_user_list.user[x])

        text_to = QLabel('Operator:')
        self.testOperator = QComboBox(self)
        for x in range(len(granta_user_list)):
            self.testOperator.addItem(granta_user_list.user[x])

        text_tl = QLabel('Test Lab:')
        self.testLab = QComboBox(self)
        for x in range(len(granta_test_lab)):
            self.testLab.addItem(granta_test_lab.lab[x])
        
        text_tg = QLabel('Test Group:')
        self.testGroup = QComboBox(self)
        for x in range(len(granta_test_group)):
            self.testGroup.addItem(granta_test_group.group[x])
        
        self.testGroup.setCurrentIndex(0)  # set default to first option

        text_loc = QLabel('Test Location:')
        self.testLoc = QComboBox(self)
        for x in range(len(granta_test_location)):
            self.testLoc.addItem(granta_test_location.location[x])
        
        text_so = QLabel('Specimen Orientation:')
        self.specimenOrient = QComboBox(self)
        for x in range(len(granta_specimen_orientation)):
            self.specimenOrient.addItem(granta_specimen_orientation.orientation[x])
    
        #self.specimenOrient.addItem('45\xb0 (TT/IP-RD)') # NOTE: when we specified discrete options inline instead of pulling from txt files, \xb0 was needed instead of the degree sign to get things to work with encodings
        self.specimenOrient.setCurrentIndex(5)   # set default to TT

        # Values are saved as the new defaults when the program closes. Try to load those, or use the standard defaults.
        try:
            self.dataSensitivity.setCurrentIndex(self.settings.value('default_dataSensitivity'))
            self.dataOriginator.setCurrentIndex(self.settings.value('default_dataOriginator'))
            self.testContact.setCurrentIndex(self.settings.value('default_testContact'))
            self.testOperator.setCurrentIndex(self.settings.value('default_testOperator'))
            self.testLab.setCurrentIndex(self.settings.value('default_testLab'))
            self.testGroup.setCurrentIndex(self.settings.value('default_testGroup'))
            self.testLoc.setCurrentIndex(self.settings.value('default_testLoc'))
            self.specimenOrient.setCurrentIndex(self.settings.value('default_specimenOrient'))
        except:
            pass
        #
        layout = QGridLayout()
        layout.addWidget(text_ds, 1, 0)
        layout.addWidget(self.dataSensitivity, 1, 1)
        layout.addWidget(text_do, 2, 0)
        layout.addWidget(self.dataOriginator, 2, 1)
        layout.addWidget(text_tc, 3, 0)
        layout.addWidget(self.testContact, 3, 1)
        layout.addWidget(text_to, 4, 0)
        layout.addWidget(self.testOperator, 4, 1)
        layout.addWidget(text_tl, 5, 0)
        layout.addWidget(self.testLab, 5, 1)
        layout.addWidget(text_tg, 6, 0)
        layout.addWidget(self.testGroup, 6, 1)
        layout.addWidget(text_loc, 7, 0)
        layout.addWidget(self.testLoc, 7, 1)
        layout.addWidget(text_so, 8, 0)
        layout.addWidget(self.specimenOrient, 8, 1)

        self.topGroupBox.setLayout(layout)

    def inputErrorDialog(self):
        msgBox = QMessageBox()
        msgBox.setIcon(QMessageBox.Information)
        msgBox.setText("Unspecified Error! Abandon Hope!")
        msgBox.setWindowTitle("Warning")
        msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        #
        returnValue = msgBox.exec()
        if returnValue == QMessageBox.Ok:
            return
        elif returnValue == QMessageBox.Cancel:
            return
#
    def cancel(self):
        print('Input aborted by user.')
        super(InputGrantaData, self).reject()
        self.close()
        sys.exit()
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
        self.OKCancel.setLayout(layout)
#
    def OKcontinue(self, event):
        self.settings.setValue('default_dataSensitivity', self.dataSensitivity.currentIndex())
        self.settings.setValue('default_dataOriginator', self.dataOriginator.currentIndex())
        self.settings.setValue('default_testContact', self.testContact.currentIndex())
        self.settings.setValue('default_testOperator', self.testOperator.currentIndex())
        self.settings.setValue('default_testLab', self.testLab.currentIndex())
        self.settings.setValue('default_testGroup', self.testGroup.currentIndex())
        self.settings.setValue('default_testLoc', self.testLoc.currentIndex())
        self.settings.setValue('default_specimenOrient', self.specimenOrient.currentIndex())
        print('GRANTA default settings saved to: ', self.settings.fileName())
        super(InputGrantaData, self).accept()

    def get_output(self):
        return (str(self.dataSensitivity.currentText()), str(self.dataOriginator.currentText()), str(self.testContact.currentText()),
            str(self.testOperator.currentText()), str(self.testLab.currentText()), str(self.testGroup.currentText()),
            str(self.testLoc.currentText()), str(self.specimenOrient.currentText()))
#
if __name__ == '__main__':
    import os, sys
    app = QApplication(sys.argv)
    main_window = InputGrantaData()
    main_window.show()
    sys.exit(app.exec_())
