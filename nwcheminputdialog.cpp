/**********************************************************************
  NWChemInputDialog - Dialog for generating NWChem input decks

  Copyright (C) 2008-2009 Marcus D. Hanwell
  Copyright (C) 2009 David C. Lonie

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.cc/>

  Avogadro is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Avogadro is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/

#include "nwcheminputdialog.h"

#include <avogadro/molecule.h>
#include <avogadro/atom.h>

#include <openbabel/mol.h>

#include <QString>
//#include <QTextStream>
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>
#include <QFile>

using namespace OpenBabel;

namespace Avogadro
{
  NWChemInputDialog::NWChemInputDialog(QWidget *parent, Qt::WindowFlags f)
    : InputDialog(parent, f), m_calculationType(OPT),
    m_theoryType(B3LYP), m_basisType(B321g),m_basisType2(B631g),m_basisType3(ccpvtz),m_tt2(PBE0),m_tt3(PBE0),
    m_output(), m_coordType(CARTESIAN), m_dirty(false), m_warned(false),nmaxiter(50),nmaxiter2(50),nmaxiter3(50),
    ntddftiter(1000),nroots(5),m_plotspin(total),m_openShell(false),m_ffreq(1.0),m_fcenter(5.0),m_fwidth(5.0),rtVis(false),
    m_rt(false),m_tmax(25.),m_dt(1.),m_rtRestart(false),m_cis(false),m_visRef(false),m_vend(25.0),m_vstart(0.),
    m_calculationType2(OPT),m_calculationType3(OPT),m_nmaxitergeom(50),m_direct(true),m_semidirect(false),m_noio(false),
    m_restartMain(false),m_dplotdens(false),m_ecp(false),m_job("molecule"),m_propall(false),m_propnbo(false),m_propesp(false),
    m_propefield(false),m_propaim(false),m_propdipole(false),m_propquad(false),m_propoct(false),m_propdens(false),
    m_propegrad(false),m_propraman(false),m_prop(false),m_dft(true),m_scf(false),m_mp2(false),m_ccsd(false)
  {
    ui.setupUi(this);
    proc = new QProcess(this);

    // Connect the GUI elements to the correct slots
    connect(ui.titleLine, SIGNAL(editingFinished()),
        this, SLOT(setTitle()));
    connect(ui.lineEdit_jobname, SIGNAL(textEdited(QString)),
        this, SLOT(setJob(QString)));
    connect(ui.calculationCombo, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setCalculation(int)));
    connect(ui.comboBox_calc2, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setCalculation2(int)));
    connect(ui.comboBox_calc3, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setCalculation3(int)));
    connect(ui.theoryCombo, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setTheory(int)));
    connect(ui.comboBox_calc2theory, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setTheory2(int)));
    connect(ui.comboBox_calc3theory, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setTheory3(int)));
    connect(ui.basisCombo, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setBasis(int)));
    connect(ui.comboBox_basis2, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setBasis2(int)));
    connect(ui.comboBox_basis3, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setBasis3(int)));
    connect(ui.comboBox_spin, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setSpin(int)));
    connect(ui.comboBox_rtspin, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setRtSpin(int)));
    connect(ui.comboBox_field, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setField(int)));
    connect(ui.comboBox_polarization, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setPolarization(int)));
    connect(ui.multiplicitySpin, SIGNAL(valueChanged(int)),
        this, SLOT(setMultiplicity(int)));
    connect(ui.chargeSpin, SIGNAL(valueChanged(int)),
        this, SLOT(setCharge(int)));
    connect(ui.coordCombo, SIGNAL(currentIndexChanged(int)),
        this, SLOT(setCoords(int)));
    connect(ui.previewText, SIGNAL(cursorPositionChanged()),
        this, SLOT(previewEdited()));
    connect(ui.generateButton, SIGNAL(clicked()),
        this, SLOT(generateClicked()));
    connect(ui.resetButton, SIGNAL(clicked()),
        this, SLOT(resetClicked()));
    connect(ui.moreButton, SIGNAL(clicked()),
        this, SLOT(moreClicked()));
    connect(ui.enableFormButton, SIGNAL(clicked()),
        this, SLOT(enableFormClicked()));
    connect(ui.checkBox_geomxyz, SIGNAL(toggled(bool)),
        this, SLOT(geomSourceChanged(bool)));
    connect(ui.checkBox_geomEditor, SIGNAL(toggled(bool)),
        this, SLOT(geomSourceChanged1(bool)));
    connect(ui.checkBox_opt2, SIGNAL(stateChanged(int)),
        this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_opt3, SIGNAL(stateChanged(int)),
        this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_xyz, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_xyz2, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_xyz3, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.spinBox_iter, SIGNAL(valueChanged(int)),
         this, SLOT(maxiterChanged(int)));
    connect(ui.spinBox_iter2, SIGNAL(valueChanged(int)),
         this, SLOT(maxiter2Changed(int)));
    connect(ui.spinBox_iter3, SIGNAL(valueChanged(int)),
         this, SLOT(maxiter3Changed(int)));
    connect(ui.checkBox_tddft, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.spinBox_tddftiter, SIGNAL(valueChanged(int)),
         this, SLOT(tddftiterChanged(int)));
    connect(ui.spinBox_nroots, SIGNAL(valueChanged(int)),
         this, SLOT(nrootsChanged(int)));
    connect(ui.spinBox_dplotnstates, SIGNAL(valueChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_tddftsinglet, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_tddfttriplet, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_mo, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.checkBox_transden, SIGNAL(stateChanged(int)),
         this, SLOT(opt2Changed(int)));
    connect(ui.lineEdit_geom, SIGNAL(editingFinished()),
          this, SLOT(setGeomFile()));
    connect(ui.checkBox_openShell, SIGNAL(toggled(bool)),
          this, SLOT(setOpenShell(bool)));
    connect(ui.checkBox_rttddft, SIGNAL(toggled(bool)),
          this, SLOT(setRttddft(bool)));
    connect(ui.checkBox_rtVis, SIGNAL(toggled(bool)),
          this, SLOT(setrtVis(bool)));
    connect(ui.doubleSpinBox_freq, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_center, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_width, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_max, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_tmax,SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_dt,SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.checkBox_restart,SIGNAL(toggled(bool)),
          this, SLOT(setRestart(bool)));
    connect(ui.checkBox_cis, SIGNAL(toggled(bool)),
          this, SLOT(setCis(bool)));
    connect(ui.checkBox_visRef, SIGNAL(toggled(bool)),
          this, SLOT(setVisRef(bool)));
    connect(ui.doubleSpinBox_refVis, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_startVis, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.doubleSpinBox_endVis, SIGNAL(valueChanged(double)),
          this, SLOT(opt2Changed(double)));
    connect(ui.spinBox_geomIter, SIGNAL(valueChanged(int)),
          this, SLOT(maxitergeomChanged(int)));
    connect(ui.checkBox_direct, SIGNAL(toggled(bool)),
          this, SLOT(setDirect(bool)));
    connect(ui.checkBox_semidirect, SIGNAL(toggled(bool)),
          this, SLOT(setSemidirect(bool)));
    connect(ui.checkBox_noio, SIGNAL(toggled(bool)),
          this, SLOT(setNoio(bool)));
    connect(ui.checkBox_restartMain, SIGNAL(toggled(bool)),
          this, SLOT(setRestartMain(bool)));
    connect(ui.checkBox_density, SIGNAL(toggled(bool)),
           this, SLOT(setDplotdens(bool)));
    connect(ui.checkBox_ecp, SIGNAL(toggled(bool)),
           this, SLOT(setEcp(bool)));
    connect(ui.lineEdit_ecp, SIGNAL(textEdited(QString)),
           this, SLOT(opt2Changed(QString)));
    connect(ui.comboBox_ecp, SIGNAL(currentIndexChanged(int)),
           this, SLOT(opt2Changed(int)));
    connect(ui.spinBox_cubemm1_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.spinBox_cubemm2_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.spinBox_cubemm3_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.spinBox_cubep1_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.spinBox_cubep2_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.spinBox_cubep3_2 ,SIGNAL(valueChanged(int)),
           this, SLOT(cube2Changed(int)));
    connect(ui.lineEdit_output, SIGNAL(textEdited(QString)),
            this, SLOT(setJob(QString)));
    connect(ui.pushButton_nwchem, SIGNAL(released()),
            this, SLOT(runNwchem()));
    connect(proc, SIGNAL(stateChanged(QProcess::ProcessState)),
            this, SLOT(jobDone(QProcess::ProcessState)));
    connect(ui.pushButton_killjob, SIGNAL(released()),
            this, SLOT(killJob()));
    connect(ui.checkBox_propertyAll ,SIGNAL(toggled(bool)),
               this, SLOT(setPropAll (bool)));
    connect(ui.checkBox_propertyEfield ,SIGNAL(toggled(bool)),
               this, SLOT(setPropEfield (bool)));
    connect(ui.checkBox_propertyAim ,SIGNAL(toggled(bool)),
               this, SLOT(setPropAim(bool)));
    connect(ui.checkBox_propertyEsp ,SIGNAL(toggled(bool)),
               this, SLOT(setPropEsp (bool)));
    connect(ui.checkBox_propertyEgrad ,SIGNAL(toggled(bool)),
               this, SLOT(setPropEgrad (bool)));
    connect(ui.checkBox_propertyDipole ,SIGNAL(toggled(bool)),
               this, SLOT(setPropDip (bool)));
    connect(ui.checkBox_propertyQuad ,SIGNAL(toggled(bool)),
               this, SLOT(setPropQuad (bool)));
    connect(ui.checkBox_propertyOct ,SIGNAL(toggled(bool)),
               this, SLOT(setPropOct (bool)));
    connect(ui.checkBox_propertyNbo ,SIGNAL(toggled(bool)),
               this, SLOT(setPropNbo (bool)));
    connect(ui.checkBox_propertyRaman ,SIGNAL(toggled(bool)),
               this, SLOT(setPropRaman (bool)));
    connect(ui.checkBox_propertyDens ,SIGNAL(toggled(bool)),
               this, SLOT(setPropDens (bool)));
    connect(ui.checkBox_prop ,SIGNAL(toggled(bool)),
               this, SLOT(setProp (bool)));

//    connect(ui.checkBox ,SIGNAL(),
//           this, SLOT(opt2Changed(int)));

    QSettings settings;
    readSettings(settings);

    // Generate an initial preview of the input deck
    setJob();
    updatePreviewText();
  }

    void NWChemInputDialog::runNwchem() {
      QString str2 = saveInputFile(ui.previewText->toPlainText(), tr("NWChem Input Deck"), QString("nw"));
      if (str2 != "") {
        QString inpFile = str2;
        QString workingDir = str2.remove(m_job+".nw");
        QString str = "mpirun -np "+QString::number(ui.spinBox_nproc->value())+" nwchem "+inpFile;

        proc->setWorkingDirectory(workingDir);
        qDebug()<<str2;
        qDebug()<<proc->workingDirectory();
        proc->setStandardOutputFile(workingDir+m_jobout);

        ui.pushButton_killjob->setEnabled(true);
        ui.pushButton_nwchem->setEnabled(false);

        proc->start(str.toLatin1());
        QString str3 = proc->readAllStandardOutput();
        QFile qFile(workingDir+m_jobout);
          if (qFile.open(QIODevice::WriteOnly)) {
            QTextStream out(&qFile); out << str3;
            qFile.close();
          }
          if ( proc->state() == QProcess::NotRunning ) {
            ui.pushButton_killjob->setEnabled(false);
            ui.pushButton_nwchem->setEnabled(true);
          };
       }
    }

    void NWChemInputDialog::jobDone(QProcess::ProcessState p) {
      if (p == 0) {
        ui.pushButton_nwchem->setEnabled(true);
        ui.pushButton_killjob->setEnabled(false);
      } else if (p == 1) {
        qDebug()<<"process running";
      }
    }

    void NWChemInputDialog::killJob() {
      QMessageBox msgBox;

      msgBox.setWindowTitle(tr("NWChem Kill Job Warning"));
      msgBox.setText(tr("Would you like to stop the current job? All progress will be saved in "+m_jobout.toLatin1()));
      msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);

      switch (msgBox.exec()) {
        case QMessageBox::Yes:
          proc->close();
          break;
        case QMessageBox::No:
          break;
        default:
          // should never be reached
          break;
      }

    }


    void NWChemInputDialog::cube2Changed(int) {
      m_dplotx = ui.spinBox_cubemm1_2->value();
      m_dploty= ui.spinBox_cubemm2_2->value();
      m_dplotz= ui.spinBox_cubemm3_2->value();
      m_dplotpx= ui.spinBox_cubep1_2->value();
      m_dplotpy= ui.spinBox_cubep2_2->value();
      m_dplotpz= ui.spinBox_cubep3_2->value();
      updatePreviewText();
    }

    void NWChemInputDialog::setEcp(bool n) {
      m_ecp = n;
      if (n) {
        ui.comboBox_ecp->setEnabled(true);
        ui.lineEdit_ecp->setEnabled(true);
      } else {
        ui.comboBox_ecp->setEnabled(false);
        ui.lineEdit_ecp->setEnabled(false);
      }
      updatePreviewText();
    }

  void NWChemInputDialog::setRestartMain(bool n) {
    m_restartMain = n;
    updatePreviewText();
  }

  void NWChemInputDialog::setDirect(bool n) {
      m_direct = n;
 //     if (m_direct)
 //       ui.checkBox_semidirect->setChecked(false);
      updatePreviewText();
  }
  void NWChemInputDialog::setSemidirect(bool n) {
      m_semidirect = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setNoio(bool n) {
      m_noio = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setProp(bool n) {
      m_prop = n;
      if (n) {
        ui.checkBox_propertyAll->setEnabled(true);
        ui.checkBox_propertyEsp->setEnabled(true);
        ui.checkBox_propertyEfield->setEnabled(true);
        ui.checkBox_propertyEgrad->setEnabled(true);
        ui.checkBox_propertyNbo->setEnabled(true);
        ui.checkBox_propertyDipole->setEnabled(true);
        ui.checkBox_propertyQuad->setEnabled(true);
        ui.checkBox_propertyOct->setEnabled(true);
        ui.checkBox_propertyRaman->setEnabled(true);
        ui.checkBox_propertyAim->setEnabled(true);
        ui.checkBox_propertyDens->setEnabled(true);
      } else {
        ui.checkBox_propertyAll->setEnabled(false);
        ui.checkBox_propertyEsp->setEnabled(false);
        ui.checkBox_propertyEfield->setEnabled(false);
        ui.checkBox_propertyEgrad->setEnabled(false);
        ui.checkBox_propertyNbo->setEnabled(false);
        ui.checkBox_propertyDipole->setEnabled(false);
        ui.checkBox_propertyQuad->setEnabled(false);
        ui.checkBox_propertyOct->setEnabled(false);
        ui.checkBox_propertyRaman->setEnabled(false);
        ui.checkBox_propertyAim->setEnabled(false);
        ui.checkBox_propertyDens->setEnabled(false);
      }
      updatePreviewText();
  }

  void NWChemInputDialog::setPropAll(bool n) {
      m_propall = n;

      updatePreviewText();
  }
  void NWChemInputDialog::setPropRaman(bool n) {
      m_propraman = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropAim(bool n) {
      m_propaim = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropEsp(bool n) {
      m_propesp = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropEfield(bool n) {
      m_propefield = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropEgrad(bool n) {
      m_propegrad = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropDip(bool n) {
      m_propdipole = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropQuad(bool n) {
      m_propquad = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropOct(bool n) {
      m_propoct = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropNbo(bool n) {
      m_propnbo = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPropDens(bool n) {
      m_propdens = n;
      updatePreviewText();
  }

  void NWChemInputDialog::setVisRef(bool n) {
    m_visRef = n;
    updatePreviewText();
  }
  void NWChemInputDialog::setDplotdens(bool n) {
    m_dplotdens = n;
    updatePreviewText();
  }

  void NWChemInputDialog::setCis(bool n) {
      m_cis = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setrtVis(bool n) {
      rtVis = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setRestart(bool n) {
      m_rtRestart = n;
      updatePreviewText();
  }


  void NWChemInputDialog::setRttddft(bool n) {
      m_rt = n;
      updatePreviewText();
  }
  void NWChemInputDialog::setField(int n) {
      m_field = (NWChemInputDialog::fieldType) n;
      updatePreviewText();
  }
  void NWChemInputDialog::setPolarization(int n) {
      m_polarization = (NWChemInputDialog::polarizationType) n;
      updatePreviewText();
  }
  void NWChemInputDialog::setOpenShell(bool n) {
      m_openShell = n;
      updatePreviewText();
  }
  void NWChemInputDialog::maxitergeomChanged(int str) {
      m_nmaxitergeom = str;
      updatePreviewText();
  }
  void NWChemInputDialog::maxiterChanged(int str) {
      nmaxiter = str;
      updatePreviewText();
  }
  void NWChemInputDialog::maxiter2Changed(int str) {
     nmaxiter2 = str;
     updatePreviewText();
  }
  void NWChemInputDialog::maxiter3Changed(int str) {
     nmaxiter3 = str;
     updatePreviewText();
  }
  void NWChemInputDialog::tddftiterChanged(int str) {
     ntddftiter = str;
     updatePreviewText();
  }
  void NWChemInputDialog::nrootsChanged(int str) {
     nroots = str;
     updatePreviewText();
  }
  void NWChemInputDialog::opt2Changed(int n)
  {
      updatePreviewText();
  }
  void NWChemInputDialog::opt2Changed(QString n)
  {
      updatePreviewText();
  }
  void NWChemInputDialog::opt2Changed(double n)
  {
      updatePreviewText();
  }

  NWChemInputDialog::~NWChemInputDialog()
  {
      QSettings settings;
      writeSettings(settings);
  }

  void NWChemInputDialog::setMolecule(Molecule *molecule)  {
    // Disconnect the old molecule first...
    if (m_molecule)
      disconnect(m_molecule, 0, this, 0);

    m_molecule = molecule;

    // Set multiplicity to the OB value
    OpenBabel::OBMol obmol = m_molecule->OBMol();
    setMultiplicity(obmol.GetTotalSpinMultiplicity());

    // Update the preview text whenever primitives are changed
    connect(m_molecule, SIGNAL(atomRemoved(Atom *)),
            this, SLOT(updatePreviewText()));
    connect(m_molecule, SIGNAL(atomAdded(Atom *)),
            this, SLOT(updatePreviewText()));
    connect(m_molecule, SIGNAL(atomUpdated(Atom *)),
            this, SLOT(updatePreviewText()));
    // Add atom coordinates
    updatePreviewText();
  }

  void NWChemInputDialog::showEvent(QShowEvent *)
  {
    updatePreviewText();
  }

  void NWChemInputDialog::updatePreviewText()
  {
    if (!isVisible())
      return;
    // Generate the input deck and display it
    if (m_dirty && !m_warned) {
      m_warned = true;
      QMessageBox msgBox;

      msgBox.setWindowTitle(tr("NWChem Input Deck Generator Warning"));
      msgBox.setText(tr("Would you like to update the preview text, losing all changes made in the NWChem input deck preview pane?"));
      msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);

      switch (msgBox.exec()) {
        case QMessageBox::Yes:
          // yes was clicked
          deckDirty(false);
          ui.previewText->setText(generateInputDeck());
          ui.previewText->document()->setModified(false);
          m_warned = false;
          break;
        case QMessageBox::No:
          // no was clicked
          m_warned = false;
          break;
        default:
          // should never be reached
          break;
      }
    }
    else if (!m_dirty) {
      ui.previewText->setText(generateInputDeck());
      ui.previewText->document()->setModified(false);
    }
  }

  void NWChemInputDialog::resetClicked()
  {
    // Reset the form to defaults
    deckDirty(false);
    ui.calculationCombo->setCurrentIndex(1);
    ui.theoryCombo->setCurrentIndex(3);
    ui.basisCombo->setCurrentIndex(2);
    ui.multiplicitySpin->setValue(0);
    ui.chargeSpin->setValue(0);
    ui.previewText->setText(generateInputDeck());
    ui.previewText->document()->setModified(false);
  }

  void NWChemInputDialog::generateClicked()
  {
    saveInputFile(ui.previewText->toPlainText(), tr("NWChem Input Deck"), QString("nw"));
  }

  void NWChemInputDialog::moreClicked()
  {
    // If the more button is clicked hide/show the preview text
    if (ui.previewText->isVisible()) {
      ui.previewText->hide();
      ui.moreButton->setText(tr("Show Preview"));
    }
    else {
      ui.previewText->show();
      ui.moreButton->setText(tr("Hide Preview"));
    }
  }

  void NWChemInputDialog::enableFormClicked()
  {
    updatePreviewText();
  }

  void NWChemInputDialog::previewEdited()
  {
    // Determine if the preview text has changed from the form generated
    if(ui.previewText->document()->isModified())
      deckDirty(true);
  }

  void NWChemInputDialog::geomSourceChanged(bool s) {
      if (ui.checkBox_geomEditor->isChecked() && ui.checkBox_geomxyz->isChecked()) {
          ui.checkBox_geomEditor->setChecked(false);
      }
      updatePreviewText();
  }
  void NWChemInputDialog::geomSourceChanged1(bool s) {
      if (ui.checkBox_geomxyz->isChecked() && ui.checkBox_geomEditor->isChecked()) {
          ui.checkBox_geomxyz->setChecked(false);
      }
      updatePreviewText();
  }


  void NWChemInputDialog::setTitle()
  {
    m_title = ui.titleLine->text();
    updatePreviewText();
  }


  void NWChemInputDialog::setJob(QString str)
  {
    m_job = ui.lineEdit_jobname->text();
    m_job = m_job.replace(" ","_");

    m_jobout = str+".out";
    ui.lineEdit_output->setText(m_jobout);

    m_molecule->setFileName(m_jobout);

    updatePreviewText();
  }

  void NWChemInputDialog::setJob()
  {
    m_job = ui.lineEdit_jobname->text();
    m_job = m_job.replace(" ","_");
    m_jobout = ui.lineEdit_jobname->text()+".out";
    ui.lineEdit_output->setText(m_jobout);

    updatePreviewText();
  }

  void NWChemInputDialog::setGeomFile() {
      m_geomFileName = ui.lineEdit_geom->text();
      updatePreviewText();
  }

  void NWChemInputDialog::setCalculation(int n)
  {
    m_calculationType = (NWChemInputDialog::calculationType) n;
    updatePreviewText();
  }
  void NWChemInputDialog::setCalculation2(int n)
  {
    m_calculationType2 = (NWChemInputDialog::calculationType) n;
    updatePreviewText();
  }
  void NWChemInputDialog::setCalculation3(int n)
  {
    m_calculationType3 = (NWChemInputDialog::calculationType) n;
    updatePreviewText();
  }

  void NWChemInputDialog::setTheory(int n)
  {
    m_theoryType = (NWChemInputDialog::theoryType) n;
    ui.basisCombo->setEnabled(true);

    if (m_theoryType == B3LYP) {
      ui.multiplicitySpin->setEnabled(true);
    } else {
      ui.multiplicitySpin->setEnabled(false);
    }

    updatePreviewText();
  }
  void NWChemInputDialog::setTheory2(int n)
  {
    m_tt2 = (NWChemInputDialog::theoryType) n;
    updatePreviewText();
  }
  void NWChemInputDialog::setTheory3(int n)
  {
    m_tt3 = (NWChemInputDialog::theoryType) n;
    updatePreviewText();
  }

  void NWChemInputDialog::setSpin(int n) {
      m_plotspin = (NWChemInputDialog::spinType) n;
      updatePreviewText();
  }
  void NWChemInputDialog::setRtSpin(int n) {
      m_plotrtspin = (NWChemInputDialog::spinType) n;
      updatePreviewText();
  }

  void NWChemInputDialog::setBasis(int n)
  {
    m_basisType = (NWChemInputDialog::basisType) n;
    updatePreviewText();
  }
  void NWChemInputDialog::setBasis2(int n)
  {
    m_basisType2 = (NWChemInputDialog::basisType) n;
    updatePreviewText();
  }
  void NWChemInputDialog::setBasis3(int n)
  {
    m_basisType3 = (NWChemInputDialog::basisType) n;
    updatePreviewText();
  }

  void NWChemInputDialog::setMultiplicity(int n)
  {
    m_multiplicity = n;
    if (ui.multiplicitySpin->value() != n) {
      ui.multiplicitySpin->setValue(n);
    }
    updatePreviewText();
  }

  void NWChemInputDialog::setOptIter(int n) {
      optiter = n;
  }
  void NWChemInputDialog::setOptIter2(int n) {
      optiter2 = n;
  }
  void NWChemInputDialog::setOptIter3(int n) {
      optiter3 = n;
  }

  void NWChemInputDialog::setCharge(int n)
  {
    m_charge = n;
    updatePreviewText();
  }

  void NWChemInputDialog::setCoords(int n)
  {
    m_coordType = (NWChemInputDialog::coordType) n;
    updatePreviewText();
  }
  QString NWChemInputDialog::getOpenShell(bool n) {
      if (n)
          return "  odft\n";
      else
          return "";
  }


  QString NWChemInputDialog::printBasisDeck(basisType t) {
      QString str;
      // Basis set
      str += "basis";

      // Need spherical keyword if using Dunning correlation consistent basis sets
      if ( t == ccpvdz || t == ccpvtz || t == augccpvdz || t == augccpvtz )
        str += " spherical";

      str += "\n";

      str += "  * library " + getBasisType(t) + '\n';
      if (m_ecp) {
        //get list of ecp atoms
        //int necp = getNEcpAtoms(ui.lineEdit_ecp->text());
        ecpType t1 = (NWChemInputDialog::ecpType) ui.comboBox_ecp->currentIndex();
        QStringList strl = ui.lineEdit_ecp->text().split(",",QString::SkipEmptyParts);
        for (int i=0; i<strl.count(); i++) {
          str += "  "+strl.at(i)+" library "+getEcpType(t1) +"\n";
        }
      }
      str += "end\n\n";

      if (m_ecp) {
        ecpType t1 = (NWChemInputDialog::ecpType) ui.comboBox_ecp->currentIndex();
        str += "ecp\n";
        QStringList strl = ui.lineEdit_ecp->text().split(",",QString::SkipEmptyParts);
        for (int i=0; i<strl.count(); i++) {
          str += "  "+strl.at(i)+" library "+getEcpType(t1) +"\n";
        }
        str += "end\n\n";
      }
      return str;
  }



  QStringList NWChemInputDialog::getEcpAtoms() {
    return ui.lineEdit_ecp->text().split(",",QString::SkipEmptyParts);
  }

  int NWChemInputDialog::getNEcpAtoms(QString str) {
    QStringList strl = str.split(",",QString::SkipEmptyParts);
    int n = strl.count();
    return n;
  }

  QString NWChemInputDialog::getEcpType(ecpType t) {
    QString str;
    switch (t) {
      case LANL2DZ_ECP:
        str += "lanl2dz_ecp";
        break;
      case CRENBL:
        str += "crenbl_ecp";
      break;
      case CRENBS:
        str += "crenbs_ecp";
      break;
      case STUTTGARTRLC:
        str += "stuttgart_rlc_ecp";
      break;
      case STUTTGARTRSC:
        str += "stuttgart_rsc_ecp";
      break;
      case SBKJCVDZ:
        str += "sbkjc_vdz_ecp";
      break;
      default:
        str += "lanl2dz_ecp";
    }
    return str;
  }

  QString NWChemInputDialog::getIo() {
    QString str;
    if (m_direct)
      str += "  direct\n";
    if (m_noio)
      str += "  noio\n";
    if (m_semidirect)
      str += "  semidirect\n";

    return str;
  }

  QString NWChemInputDialog::printTheory(theoryType t ) {
      QString str;
      switch (t)
        {
        case PBE0:
         str += "dft\n";
         str += getIo();
         str += getOpenShell(m_openShell)+"  xc pbe0\n  mulliken\n  maxiter "+QString::number(nmaxiter)+"\n  mult " + QString::number(m_multiplicity) + "\nend\n\n";
        break;
        case M062X:
        str += "dft\n";
        str += getIo();
        str += getOpenShell(m_openShell)+"  xc pbe0\n  mulliken\n  maxiter "+QString::number(nmaxiter)+"\n  mult " + QString::number(m_multiplicity) + "\nend\n\n";
        break;
        case B3LYP:
        str += "dft\n";
        str += getIo();
        str += getOpenShell(m_openShell)+"  xc pbe0\n  mulliken\n  maxiter "+QString::number(nmaxiter)+"\n  mult " + QString::number(m_multiplicity) + "\nend\n\n";
          break;
        case MP2:
          str += "mp2\n";
          str += "  # Exclude core electrons from MP2 treatment\n";
          str += "  freeze atomic\n";
          str += "end\n\n";
          break;
        case CCSD:
          str += "ccsd\n";
          str += "  # Exclude core electrons from CCSD treatment\n";
          str += "  freeze atomic\n";
          str += "end\n\n";
          break;
        default:
        case RHF:
            break;
        }
    return str;
  }

  QString NWChemInputDialog::printTask(theoryType t) {
    QString str = "task ";
    switch (t) {
      case PBE0:
        str += "dft ";
        break;
      case M062X:
        str += "dft ";
        break;
      case B3LYP:
        str += "dft ";
        break;
      case CCSD:
        str += "ccsd ";
        break;
      case MP2:
        str += "mp2 ";
        break;
      default:
      case RHF:
        str += "scf ";
        break;
    }
    return str;
  }

  QString NWChemInputDialog::printProp() {
    QString str;
    str += "property\n";
    if (m_propaim) str += "  aimfile\n";
    if (m_propall) str += "  all\n";
    if (m_propefield ) str += "  efield\n";
    if (m_propegrad) str += "  efieldgrad\n";
    if (m_propesp) str += "  esp\n";

    str += "end \n\n";
    str += "task ";
    if (m_dft) str += "dft ";
    if (m_scf) str += "scf ";
    str += "property\n\n";

    return str;
  }

  QString NWChemInputDialog::generateInputDeck()
  {
    // Generate an input deck based on the settings of the dialog
    QString buffer;
    QTextStream mol(&buffer);
    ui.lineEdit_output->setText(m_jobout);

    //Print header
    mol << "######################################################\n";
    mol << "#   NWChem Input file created using Avogadro \n";
    mol << "#   Extended functionality provided by CTC (2019)\n";
    mol << "######################################################\n";

    // Print input in output
    mol << "echo\n";

    // Job - CTC
    if (m_rtRestart || m_restartMain) {
      mol << "restart "<< m_job <<" \n";
    } else {
      mol << "start "<< m_job <<" \n";
    }

    // Title
    mol << "title \"" << m_title << "\"\n";

    // Ecce file
    mol << "ecce_print "<< m_job <<"_ecce.out \n";

    // Now for the charge
    mol << "charge " << m_charge << "\n\n";

    // Geometry specification
    if (!m_rtRestart && !m_restartMain) {
    mol << "geometry \""<<m_job<< "\" units angstroms print";

    // Now to output the actual molecular coordinates
    // Cartesian coordinates
    if (m_molecule && m_coordType == CARTESIAN)
    {
      QTextStream mol(&buffer);
      mol << " xyz autosym\n";

      if (ui.checkBox_geomEditor->isChecked() ) {
        QList<Atom *> atoms = m_molecule->atoms();
        foreach (Atom *atom, atoms) {
            mol << qSetFieldWidth(4) << right
                << QString(OpenBabel::etab.GetSymbol(atom->atomicNumber()))
                << qSetFieldWidth(15) << qSetRealNumberPrecision(5) << forcepoint
                << fixed << right << atom->pos()->x() << atom->pos()->y()
                << atom->pos()->z()
                << qSetFieldWidth(0) << '\n';
        }
      }
    }

    if (ui.checkBox_geomxyz->isChecked() ) {
        mol << "  load "<<m_geomFileName<<"\n";
    }

    // Z-matrix
    else if (m_molecule && m_coordType == ZMATRIX)
    {
      QTextStream mol(&buffer);
      mol.setFieldAlignment(QTextStream::AlignAccountingStyle);
      mol << "\n zmatrix\n";
      OBAtom *a, *b, *c;
      double r, w, t;

      /* Taken from OpenBabel's gzmat file format converter */
      std::vector<OBInternalCoord*> vic;
      vic.push_back((OpenBabel::OBInternalCoord*)NULL);
      OpenBabel::OBMol obmol = m_molecule->OBMol();
      FOR_ATOMS_OF_MOL(atom, &obmol)
        vic.push_back(new OpenBabel::OBInternalCoord);
      CartesianToInternal(vic, obmol);

      FOR_ATOMS_OF_MOL(atom, &obmol)
      {
        a = vic[atom->GetIdx()]->_a;
        b = vic[atom->GetIdx()]->_b;
        c = vic[atom->GetIdx()]->_c;

        mol << qSetFieldWidth(3) << QString(etab.GetSymbol(atom->GetAtomicNum()));

        if (atom->GetIdx() > 1)
          mol << qSetFieldWidth(0) << "  " << qSetFieldWidth(3) << QString::number(a->GetIdx())
              << qSetFieldWidth(0) << "  "<< qSetFieldWidth(4) << QString("r") + QString::number(atom->GetIdx());

        if (atom->GetIdx() > 2)
          mol << qSetFieldWidth(0) << "  " << qSetFieldWidth(3) << QString::number(b->GetIdx())
              << qSetFieldWidth(0) << "  "<< qSetFieldWidth(4) << QString("a") + QString::number(atom->GetIdx());

        if (atom->GetIdx() > 3)
          mol << qSetFieldWidth(0) << "  " << qSetFieldWidth(3) << QString::number(c->GetIdx())
              << qSetFieldWidth(0) << "  "<< qSetFieldWidth(4) << QString("d") + QString::number(atom->GetIdx());

        mol << qSetFieldWidth(0) << '\n';
      }

      mol << " variables\n";
      FOR_ATOMS_OF_MOL(atom, &obmol)
      {
        r = vic[atom->GetIdx()]->_dst;
        w = vic[atom->GetIdx()]->_ang;
        if (w < 0.0)
          w += 360.0;
        t = vic[atom->GetIdx()]->_tor;
        if (t < 0.0)
          t += 360.0;
        if (atom->GetIdx() > 1)
          mol << "   r" << atom->GetIdx() << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right
              << r << qSetFieldWidth(0) << '\n';
        if (atom->GetIdx() > 2)
          mol << "   a" << atom->GetIdx() << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right
              << w << qSetFieldWidth(0) << '\n';
        if (atom->GetIdx() > 3)
          mol << "   d" << atom->GetIdx() << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right
              << t << qSetFieldWidth(0) << '\n';
      }
      mol << " end\n";
      foreach (OpenBabel::OBInternalCoord *c, vic)
        delete c;
    }
    // Compact ZMatrix
    else if (m_molecule && m_coordType == ZMATRIX_COMPACT)
    {
      QTextStream mol(&buffer);
      mol << " zmatrix\n";
      OBAtom *a, *b, *c;
      double r, w, t;

      /* Taken from OpenBabel's gzmat file format converter */
      std::vector<OBInternalCoord*> vic;
      vic.push_back((OpenBabel::OBInternalCoord*)NULL);
      OpenBabel::OBMol obmol = m_molecule->OBMol();
      FOR_ATOMS_OF_MOL(atom, &obmol)
        vic.push_back(new OpenBabel::OBInternalCoord);
      CartesianToInternal(vic, obmol);

      FOR_ATOMS_OF_MOL(atom, &obmol)
      {
        a = vic[atom->GetIdx()]->_a;
        b = vic[atom->GetIdx()]->_b;
        c = vic[atom->GetIdx()]->_c;
        r = vic[atom->GetIdx()]->_dst;
        w = vic[atom->GetIdx()]->_ang;
        if (w < 0.0)
          w += 360.0;
        t = vic[atom->GetIdx()]->_tor;
        if (t < 0.0)
          t += 360.0;

        mol << qSetFieldWidth(4) << right
            << QString(etab.GetSymbol(atom->GetAtomicNum())
                       + QString::number(atom->GetIdx()));
        if (atom->GetIdx() > 1)
          mol << qSetFieldWidth(6) << right
              << QString(etab.GetSymbol(a->GetAtomicNum())
                         + QString::number(a->GetIdx())) << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right << r;
        if (atom->GetIdx() > 2)
          mol << qSetFieldWidth(6) << right
                 << QString(etab.GetSymbol(b->GetAtomicNum())
                         + QString::number(b->GetIdx())) << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right << w;
        if (atom->GetIdx() > 3)
          mol << qSetFieldWidth(6) << right
              << QString(etab.GetSymbol(c->GetAtomicNum())
                         + QString::number(c->GetIdx())) << qSetFieldWidth(15)
              << qSetRealNumberPrecision(5) << forcepoint << fixed << right << t;
        mol << qSetFieldWidth(0) << '\n';
      }
      foreach (OpenBabel::OBInternalCoord *c, vic)
        delete c;
    }
    mol << "end\n\n";


    //Set geometry
    mol << "set geometry \""<<m_job<<"\" \n\n";

    /*****************
     * Basis
     * ***************/
    mol << printBasisDeck(m_basisType);

    /*****************
     * Driver
     *****************/
    if (m_calculationType == OPT) {
    mol << "driver\n  maxiter  "<< m_nmaxitergeom <<"\n";
        if (ui.checkBox_xyz->isChecked() ) {
            mol << "  xyz\n";
        }
        mol<<"end\n\n" ;
    }

    /*****************
     * Theory
     * ***************/
    mol << printTheory(m_theoryType);

    /*****************
     * Task
     * ***************/
    mol << printTask(m_theoryType) ;
    mol << getCalculationType(m_calculationType) <<"\n\n"<< endl;

    } //end if not restart

    if (m_restartMain) {
      /*****************
       * Basis
       * ***************/
      mol << printBasisDeck(m_basisType);
    }

    /**************************************************************************/
    // Second Calculation
    /**************************************************************************/

    //Create DFT block
    if (ui.checkBox_opt2->isChecked() ) {
        ui.checkBox_opt3->setEnabled(true);
        mol << "###############################\n";
        if (m_calculationType2 == OPT)
          mol << "# Additional Optimization \n";
        else if (m_calculationType2 == SP)
          mol << "# Additional Single Point \n";
        else if (m_calculationType2 == FREQ)
          mol << "# Additional Frequency \n";
        mol << "###############################\n"<<endl;

        if (m_calculationType2 == OPT) {
        mol << "driver\n  maxiter  "<< nmaxiter2 <<"\n";
            if (ui.checkBox_xyz2->isChecked() ) {
                mol << "  xyz\n";
            }
            mol<<"end\n\n" ;
        }

        mol << printTheory(m_tt2);


        /********************
         * Print Basis library
         */
          // Need spherical keyword if using Dunning correlation consistent basis sets
        mol << printBasisDeck(m_basisType2);
        mol << printTask(m_tt2);
        mol << getCalculationType(m_calculationType2) <<"\n"<<endl;
    }

    /**************************************************************************/
    // Third Calculation
    /**************************************************************************/
    if (ui.checkBox_opt3->isChecked() ) {
      mol << "###############################\n";
      if (m_calculationType3 == OPT)
        mol << "# Additional Optimization \n";
      else if (m_calculationType3 == SP)
        mol << "# Additional Single Point \n";
      else if (m_calculationType3 == FREQ)
        mol << "# Additional Frequency \n";
      mol << "###############################\n"<<endl;

      if (m_calculationType3 == OPT) {
        mol << "driver\n  maxiter  "<< nmaxiter3 <<"\n";
        if (ui.checkBox_xyz3->isChecked() ) {
          mol << "  xyz\n";
        }
        mol<<"end\n\n" ;
      }
      mol << printTheory(m_tt3);

      /********************
       * Print Basis
       */
      mol << printBasisDeck(m_basisType3);
      mol << printTask(m_tt3);
      mol << getCalculationType(m_calculationType3) <<"\n"<<endl;

    }

    /***************************************
     * Properties
     * *************************************/
    if (m_prop) {
      mol << printProp() <<endl;
    }

    // Add TDDFT if needed
    if (ui.checkBox_tddft->isChecked()) {
        mol << "##########################" <<endl;
        mol << "# TDDFT" << endl;
        mol << "##########################\n" << endl;

        mol<< "tddft\n  nroots "<<nroots<<"\n  civecs\n  transden\n  maxiter "<<ntddftiter<<"\n";
        if (!ui.checkBox_tddftsinglet->isChecked()) {
            mol<<"  nosinglet\n";
        }
        if (!ui.checkBox_tddfttriplet->isChecked()) {
            mol<<"  notriplet\n";
        }
        if (m_cis)
            mol<<"  cis\n";

        mol<<"end\n\n";
        mol<<"task tddft\n\n";
    }

    // Add Dplot if needed
    if (ui.checkBox_mo->isChecked()) {
        int nmos = ui.spinBox_frontier->value();
        for (int i=0; i< 2*nmos+2; i++) {
            int moindex = ui.spinBox_homo->value()-nmos + i;
            QString cubeName = m_job+"-mo-"+QString::number(moindex)+".cube" ;
            mol << printDplotMos(cubeName, moindex);
        }
    }
    if (ui.checkBox_density->isChecked() ) {
        QString cubeName = m_job+"-density.cube" ;
        mol << printDplotDens(cubeName);
    }

    if (ui.checkBox_transden->isChecked()) {
        int ntransden = ui.spinBox_dplotnstates->value();
        for (int i=0; i<ntransden; i++) {
            QString cubeName = m_job+"-tdens-"+QString::number(i+1)+".cube";
            mol << printDplotTransden(cubeName, i+1);
        }
    }

    /***************************
     * RT-TDDFT
     * *************************/
    if (m_rt) {
        mol << printRttddft();
    }

    return buffer;
  }

/**********************************************************/
  // FUNCTIONS
/**********************************************************/

  //Function to convert between atomic and SI units
  // dir: 0 returns SI
  // dir: 1 returns au
/*
  Quantity	Conversion
  Time 1 au  = 0.02419 fs
  Length 1 au = 0.5292 A
  Energy 1 au = 27.2114 eV
  Electric field 1 au = 514.2 V/nm
  Dipole moment	1 au = 2.542 D
*/
  double NWChemInputDialog::convertUnits(int which, int dir, double val) {
      switch(which) {
        case(0): //time
          if (dir == 0)
            return val*0.02419;
          else if (dir == 1)
              return val/0.02419;
          break;
        case(1): //field
        if (dir == 0)
          return val*514200.;
        else if (dir == 1)
            return val/514200.;
        break;
        case(2): //energy
        if (dir == 0)
          return val*27.2114;
        else if (dir == 1)
            return val/27.211;
        break;
      }
  }

  QString NWChemInputDialog::printRttddft() {
      QString str;
      str += "##########################\n";
      str +="# RT-TDDFT \n";
      str += "##########################\n\n";
      str += "unset rt_tddft:* \n\n";
      str += "rt_tddft \n" ;
      str += "  tmax "+QString::number(convertUnits(0,1,ui.doubleSpinBox_tmax->value()))+"\n";
      str += "  dt "+ QString::number(convertUnits(0,1,ui.doubleSpinBox_dt->value()))+"\n";
      str += "  print *\n\n";
      str += "  field \"driver\"\n";
      str += "    type "+getfieldType(m_field)+"\n";
      str += "    polarization "+getpolarizationType(m_polarization)+"\n";
      str += "    frequency "+QString::number(convertUnits(2,1,ui.doubleSpinBox_freq->value()))+"\n";
      str += "    center "+QString::number(convertUnits(0,1,ui.doubleSpinBox_center->value()))+"\n";
      str += "    width "+QString::number(convertUnits(0,1,ui.doubleSpinBox_width->value()))+"\n";
      str += "    max "+QString::number(convertUnits(1,1,ui.doubleSpinBox_max->value()))+"\n";
      str += "  end\n\n";

      str += "  excite \""+m_job+"\" with \"driver\"\n\n";

      if (rtVis) {
          str += "  visualization\n";
          str += "    tstart "+QString::number(convertUnits(0,1,ui.doubleSpinBox_startVis->value()))+"\n";
          str += "    tend "+QString::number(convertUnits(0,1,ui.doubleSpinBox_endVis->value()))+"\n";
          if (m_visRef)
            str += "    treference "+QString::number(convertUnits(0,1,ui.doubleSpinBox_refVis->value()))+"\n";
          str += "    dplot\n";
          str += "  end\n";
      }
      str += "end\n\n";

      if (rtVis) {
          str += printDplotRt();
      }
      str += "task dft rt_tddft\n";
      return str;
  }

  QString NWChemInputDialog::printDplotRt() {
      QString str;
      str += "####################################\n";
      str +="# DPLOT for RT-TDDFT visualization\n";
      str += "####################################\n";
      str += "dplot\n" ;
      str += " limitxyz\n  -"+QString::number(ui.spinBox_cubemm1_2->value())+" "+QString::number(ui.spinBox_cubemm1_2->value()) +" "+QString::number(ui.spinBox_cubep1_2->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm2_2->value())+" "+QString::number(ui.spinBox_cubemm2_2->value())+" "+QString::number(ui.spinBox_cubep2_2->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm3_2->value())+" "+QString::number(ui.spinBox_cubemm3_2->value())+" "+QString::number(ui.spinBox_cubep3_2->value())+"\n";
      str += " gaussian\n spin "+ getSpinType(m_plotrtspin )+"\n";
      str += "end\n\n";

      return str;
  }
  QString NWChemInputDialog::printDplotMos(QString str1, int i) {
      QString str;
      str += "##########################\n";
      str +="# DPLOT\n";
      str += "##########################\n";
      str += "dplot\n" ;
      str += " limitxyz\n  -"+QString::number(ui.spinBox_cubemm1->value())+" "+QString::number(ui.spinBox_cubemm1->value()) +" "+QString::number(ui.spinBox_cubep1->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubep2->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubep3->value())+"\n";
      str += " gaussian\n spin "+ getSpinType(m_plotspin ) +"\n vectors "+m_job+".movecs\n";
      str += " orbitals view; 1; "+QString::number(i)+"\n";
      str += " output "+str1+"\nend\n";
      str += "task dplot\n\n";

      return str;
  }
  QString NWChemInputDialog::printDplotDens(QString str1) {
      QString str;
      str += "##########################\n";
      str +="# DPLOT\n";
      str += "##########################\n";
      str += "dplot\n" ;
      str += " limitxyz\n  -"+QString::number(ui.spinBox_cubemm1->value())+" "+QString::number(ui.spinBox_cubemm1->value()) +" "+QString::number(ui.spinBox_cubep1->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubep2->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubep3->value())+"\n";
      str += " gaussian\n spin "+ getSpinType(m_plotspin ) +"\n vectors "+m_job+".movecs\n";
      str += " output "+str1+"\nend\n";
      str += "task dplot\n\n";

      return str;
  }
  QString NWChemInputDialog::printDplotTransden(QString str1, int i) {
      QString str;
      str += "##########################\n";
      str +="# DPLOT\n";
      str += "##########################\n";
      str += "dplot\n" ;
      str += " limitxyz\n  -"+QString::number(ui.spinBox_cubemm1->value())+" "+QString::number(ui.spinBox_cubemm1->value()) +" "+QString::number(ui.spinBox_cubep1->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubemm2->value())+" "+QString::number(ui.spinBox_cubep2->value())+"\n";
      str += "  -"+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubemm3->value())+" "+QString::number(ui.spinBox_cubep3->value())+"\n";
      str += " gaussian\n spin "+ getSpinType(m_plotspin ) +"\n vectors "+m_job+".movecs\n";
      if (i<10)
        str += " densmat "+m_job+".tdens_00000"+QString::number(i)+"\n";
      else if (i>9 && i<100)
        str += " densmat "+m_job+".tdens_0000"+QString::number(i)+"\n";
      str += " output "+str1+"\nend\n";
      str += "task dplot\n\n";

      return str;
  }

  QString NWChemInputDialog::getfieldType(fieldType t) {
    switch (t) {
        case KICK:
            return "kick";
        case GAUSSIAN:
            return "gaussian";
        default:
            return "gaussian";
    }
  }
  QString NWChemInputDialog::getpolarizationType(polarizationType t) {
    switch (t) {
        case X:
            return "x";
        case Y:
            return "y";
        case Z:
            return "z";
        default:
            return "x";
    }
  }

  QString NWChemInputDialog::getCalculationType(calculationType t)
  {
    // Translate the enum to text for the output generation
    switch (t)
      {
      case SP:
        return "energy";
      case OPT:
        return "optimize";
      case FREQ:
        return "freq";
      default:
        return "";
      }
  }

  QString NWChemInputDialog::getTheoryType(theoryType t)
  {
    // Translate the enum to text for the output generation
    switch (t)
    {//   enum theoryType{RHF, B3LYP, B3LYP5, EDF1, M062X, MP2, CCSD}
      case RHF:
        return "RHF";
      case B3LYP:
        return "B3LYP";
      case MP2:
        return "MP2";
      case CCSD:
        return "CCSD";
      case PBE0:
        return "PBE0";
      case M062X:
        return "M062X";
      default:
        return "PBE0";
    }
  }

  QString NWChemInputDialog::getBasisType(basisType t)
  {
    // Translate the enum to text for the output generation
    switch (t)
    {
      case STO3G:
        return "STO-3g";
      case B321g:
        return "3-21g";
      case B631g:
         return "6-31g";
      case B631gp:
        return "6-31g*";
      case B631plusgp:
        return "6-31+g*";
      case B6311g:
        return "6-311g";
      case B6311gp:
        return "6-311g*";
      case ccpvdz:
        return "cc-pvdz";
      case ccpvtz:
        return "cc-pvtz";
      case augccpvdz:
        return "aug-cc-pvdz";
      case augccpvtz:
        return "aug-cc-pvtz";
      case LANL2DZ:
        return "LANL2DZ ECP";
      default:
        return "3-21g";
    }
  }

  QString NWChemInputDialog::getSpinType(spinType t) {
      switch (t) {
        case total:
          return "total";
        case alpha:
          return "alpha";
        case beta:
          return "beta";
        case spindens:
          return "spindens";
        default:
          return "total";
      }
  }

  void NWChemInputDialog::deckDirty(bool dirty)
  {
    m_dirty = dirty;
    ui.titleLine->setEnabled(!dirty);
    ui.calculationCombo->setEnabled(!dirty);
    ui.theoryCombo->setEnabled(!dirty);
    ui.basisCombo->setEnabled(!dirty);
    ui.multiplicitySpin->setEnabled(!dirty);
    ui.chargeSpin->setEnabled(!dirty);
    ui.enableFormButton->setEnabled(dirty);
  }

  void NWChemInputDialog::readSettings(QSettings& settings)
  {
    m_savePath = settings.value("nwchem/savepath").toString();
  }
  
  void NWChemInputDialog::writeSettings(QSettings& settings) const
  {
    settings.setValue("nwchem/savepath", m_savePath);
  }
}

