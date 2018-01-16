import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QLabel, QLineEdit, QPushButton,
                            QGridLayout, QGroupBox, QVBoxLayout, QListWidget,
                            QCheckBox, QMessageBox)
import matplotlib.pyplot as plt

from ptkinetic import Kinetic

class App(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('ptkinetic program interface')
        self.setFixedSize(1000,500)
        self.statusBar().showMessage('Thien-Phuc Tu-Nguyen')
        self.initUI()
        self.show()
        self.kinetic = Kinetic()
        self.chemical_name = []
        self.reactants = []
        self.products = []
        self.chemical_graph = []

    def initUI(self):
        label = QLabel(self)
        label.setText('Chemical name')
        label.move(20,20)
        label.resize(100,30)

        self.line_chemical = QLineEdit(self)
        self.line_chemical.setToolTip('Enter the chemical name here')
        self.line_chemical.move(20,60)
        self.line_chemical.resize(100,30)

        label = QLabel(self)
        label.setText('Concentration')
        label.move(20,100)
        label.resize(100,30)

        self.line_concentration = QLineEdit(self)
        self.line_concentration.move(20, 140)
        self.line_concentration.resize(100,30)

        self.check_stable = QCheckBox(self)
        self.check_stable.setText('remain stable')
        self.check_stable.move(20, 180)
        self.check_stable.resize(120,30)
        self.check_stable.setChecked(False)

        self.button_chemical = QPushButton(self)
        self.button_chemical.setText('Add chemical')
        self.button_chemical.move(20, 220)
        self.button_chemical.resize(100, 30)
        self.button_chemical.clicked.connect(self.add_chemical)

        self.list_chemical = QListWidget(self)
        self.list_chemical.move(150, 30)
        self.list_chemical.resize(100, 220)

        label = QLabel(self)
        label.setText('Rate constant')
        label.move(270, 20)
        label.resize(90, 30)

        self.line_rate = QLineEdit(self)
        self.line_rate.move(270, 70)
        self.line_rate.resize(90, 30)

        self.button_reactant = QPushButton(self)
        self.button_reactant.setText('Add reactant')
        self.button_reactant.move(270,120)
        self.button_reactant.resize(90,30)
        self.button_reactant.clicked.connect(self.add_reactant)

        self.button_product = QPushButton(self)
        self.button_product.setText('Add product')
        self.button_product.move(270,170)
        self.button_product.resize(90,30)
        self.button_product.clicked.connect(self.add_product)

        self.button_reaction = QPushButton(self)
        self.button_reaction.setText('Add reaction')
        self.button_reaction.move(270, 220)
        self.button_reaction.resize(90,30)
        self.button_reaction.clicked.connect(self.add_reaction)

        label = QLabel(self)
        label.setText('List of reactants')
        label.move(380, 30)
        label.resize(100,30)

        self.list_reactant = QListWidget(self)
        self.list_reactant.move(380, 70)
        self.list_reactant.resize(100, 180)

        label = QLabel(self)
        label.setText('List of products')
        label.move(500, 30)
        label.resize(100,30)

        self.list_product = QListWidget(self)
        self.list_product.move(500, 70)
        self.list_product.resize(100, 180)

        label = QLabel(self)
        label.setText('List of reactions')
        label.move(620, 30)
        label.resize(300,30)

        self.list_reaction = QListWidget(self)
        self.list_reaction.move(620, 70)
        self.list_reaction.resize(300, 180)

        self.button_init = QPushButton(self)
        self.button_init.setText('Initialize')
        self.button_init.move(20, 280)
        self.button_init.resize(100, 30)
        self.button_init.clicked.connect(self.init)

        label = QLabel(self)
        label.setText('Steps per cycle')
        label.move(20, 330)
        label.resize(100, 30)

        self.line_step = QLineEdit(self)
        self.line_step.move(20, 380)
        self.line_step.resize(100,30)

        label = QLabel(self)
        label.setText('Number of Cycle')
        label.move(140, 280)
        label.resize(100, 30)

        self.line_cycle = QLineEdit(self)
        self.line_cycle.move(140, 320)
        self.line_cycle.resize(100,30)
        self.line_cycle.setEnabled(False)

        label = QLabel(self)
        label.setText('Time per cycle')
        label.move(140, 360)
        label.resize(100, 30)

        self.line_time = QLineEdit(self)
        self.line_time.move(140, 410)
        self.line_time.resize(100,30)
        self.line_time.setEnabled(False)

        self.button_run = QPushButton(self)
        self.button_run.setText('Run')
        self.button_run.move(260, 320)
        self.button_run.resize(100,30)
        self.button_run.clicked.connect(self.run)
        self.button_run.setEnabled(False)

        label = QLabel(self)
        label.setText('File name')
        label.move(380, 280)
        label.resize(100,30)

        self.line_file = QLineEdit(self)
        self.line_file.move(380, 320)
        self.line_file.resize(100, 30)
        self.line_file.setEnabled(False)

        self.button_file = QPushButton(self)
        self.button_file.setText('Save')
        self.button_file.move(380, 360)
        self.button_file.resize(100,30)
        self.button_file.setEnabled(False)
        self.button_file.clicked.connect(self.save)

        self.list_graph = QListWidget(self)
        self.list_graph.move(500, 280)
        self.list_graph.resize(100, 200)
        self.list_graph.setEnabled(False)

        self.button_add = QPushButton(self)
        self.button_add.setText('Add')
        self.button_add.move(630, 280)
        self.button_add.resize(100,30)
        self.button_add.setEnabled(False)
        self.button_add.clicked.connect(self.add_graph)

        self.button_remove = QPushButton(self)
        self.button_remove.setText('Remove')
        self.button_remove.move(630, 320)
        self.button_remove.resize(100,30)
        self.button_remove.setEnabled(False)
        self.button_remove.clicked.connect(self.remove_graph)

        self.button_graph = QPushButton(self)
        self.button_graph.setText('Plot')
        self.button_graph.move(630, 360)
        self.button_graph.resize(100,30)
        self.button_graph.setEnabled(False)
        self.button_graph.clicked.connect(self.plot)

    def add_chemical(self):
        conc = None
        if self.line_chemical.text() == '':
            QMessageBox.warning(self, 'Value Error','Please enter the chemical name',
                        QMessageBox.Ok, QMessageBox.Ok)
            return
        if self.line_concentration.text() == '':
            QMessageBox.warning(self, 'Value Error','Please enter the initial concentration',
                        QMessageBox.Ok, QMessageBox.Ok)
            return
        try:
            conc = float(self.line_concentration.text())
        except ValueError:
            QMessageBox.warning(self, 'Value Error','The concentration you enter is not correct',
                        QMessageBox.Ok, QMessageBox.Ok)
            return
        try:
            self.kinetic.add_chemical(self.line_chemical.text(), float(self.line_concentration.text()),self.check_stable.isChecked())
        except ValueError:
            QMessageBox.warning(self, 'Value Error','The chemical already existed',
                        QMessageBox.Ok, QMessageBox.Ok)
            return
        self.list_chemical.addItem(self.line_chemical.text() + ' ' + self.line_concentration.text())
        self.chemical_name.append(self.line_chemical.text())
        self.line_chemical.setText('')
        self.line_concentration.setText('')
        self.check_stable.setChecked(False)

    def add_reactant(self):
        if self.list_chemical.currentRow() == -1:
            QMessageBox.warning(self, 'Value Error', 'Please select an chemical',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        ind = self.list_chemical.currentRow()
        self.list_reactant.addItem(self.chemical_name[ind])
        self.reactants.append(self.chemical_name[ind])


    def add_product(self):
        if self.list_chemical.currentRow() == -1:
            QMessageBox.warning(self, 'Value Error', 'Please select an chemical',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        ind = self.list_chemical.currentRow()
        self.list_product.addItem(self.chemical_name[ind])
        self.products.append(self.chemical_name[ind])

    def add_reaction(self):
        k = None
        if self.line_rate.text() == '':
            QMessageBox.warning(self, 'Value Error', 'Please enter the rate constant',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        if len(self.reactants) == 0:
            QMessageBox.warning(self, 'Value Error', 'Please add reactants before adding',
                                QMessageBox.Ok, QMessageBox.Ok)
        if len(self.products) == 0:
            QMessageBox.warning(self, 'Value Error', 'Please add products before adding',
                                QMessageBox.Ok, QMessageBox.Ok)
        try:
            k = float(self.line_rate.text())
        except ValueError:
            QMessageBox.warning(self, 'Value Error', 'The rate constant is not exceptable',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        self.kinetic.add_reaction(self.reactants, self.products, k)
        s = self.reactants[0]
        for i in range(1, len(self.reactants)):
            s += '+'+self.reactants[i]
        s += '->'
        s += self.products[0]
        for i in range(1,len(self.products)):
            s += '+' + self.products[i]
        s += ' ' + str(k)
        self.list_reaction.addItem(s)
        self.reactants = []
        self.products = []
        self.list_reactant.clear()
        self.list_product.clear()
        self.line_rate.setText('')

    def init(self):
        if self.list_reaction.count() == 0:
            QMessageBox.warning(self, 'Value Error', 'Please declare reactions',
                                QMessageBox.Ok, QMessageBox.Ok)
            return

        if self.line_step.text() == '':
            QMessageBox.warning(self, 'Value Error','Please enter the number of step',
                        QMessageBox.Ok, QMessageBox.Ok)
            return
        try:
            conc = int(self.line_step.text())
        except ValueError:
            QMessageBox.warning(self, 'Value Error','The number of step you enter is not correct',
                        QMessageBox.Ok, QMessageBox.Ok)
            return

        self.kinetic.init(int(self.line_step.text()))
        self.line_chemical.setEnabled(False)
        self.line_concentration.setEnabled(False)
        self.check_stable.setEnabled(False)
        self.button_chemical.setEnabled(False)
        self.button_reactant.setEnabled(False)
        self.button_product.setEnabled(False)
        self.button_reaction.setEnabled(False)
        self.list_reactant.setEnabled(False)
        self.list_product.setEnabled(False)
        self.line_rate.setEnabled(False)
        self.list_reaction.setEnabled(False)
        self.button_init.setEnabled(False)
        self.line_step.setEnabled(False)
        self.line_cycle.setEnabled(True)
        self.line_time.setEnabled(True)
        self.button_run.setEnabled(True)

    def run(self):
        cycle = None
        delta = None
        if self.line_cycle.text() == '':
            QMessageBox.warning(self, 'Value Error', 'Please enter the number of cycles',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        try:
            cycle = int(self.line_cycle.text())
        except ValueError:
            QMessageBox.warning(self, 'Value Error', 'The number of cycles is not correct',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        if self.line_time.text() == '':
            QMessageBox.warning(self, 'Value Error', 'Please enter the time step',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        try:
            delta = float(self.line_time.text())
        except ValueError:
            QMessageBox.warning(self, 'Value Error', 'The time step is not correct',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        self.kinetic.run(cycle, delta)
        self.line_file.setEnabled(True)
        self.button_file.setEnabled(True)
        self.list_graph.setEnabled(True)
        self.button_add.setEnabled(True)
        self.button_remove.setEnabled(True)
        self.button_graph.setEnabled(True)

    def save(self):
        if self.line_file.text() == '':
            QMessageBox.warning(self, 'Value Error', 'Please enter the file name',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        self.kinetic.save(self.line_file.text())

    def add_graph(self):
        index = self.list_chemical.currentRow()
        chosen = self.chemical_name[index]
        if chosen not in self.chemical_graph:
            self.list_graph.addItem(chosen)
            self.chemical_graph.append(chosen)

    def remove_graph(self):
        if self.list_graph.currentRow() == -1:
            QMessageBox.warning(self, 'Value Error', 'Please select an chemical to remove',
                                QMessageBox.Ok, QMessageBox.Ok)
            return
        index = self.list_graph.currentRow()
        self.list_graph.takeItem(index)
        self.chemical_graph.pop(index)
        return

    def plot(self):
        self.kinetic.plot(self.chemical_graph)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = App()
    sys.exit(app.exec_())
