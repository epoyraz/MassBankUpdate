import os, sys
import decimal
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1], f))]

print onlyfiles
# from collections import namedtuple
# from pybel import *
from pyopenms import *

MS_LEVEL = "MS2"

#filelist = onlyfiles
filenames = onlyfiles
#filelist.close()

mzmlfile = MzMLFile()
msexp = MSExperiment()

mb2hmdbfile = open("MB2HMDBmapping.csv", "r")
mb2hmdblines = mb2hmdbfile.readlines()
mb2hmdbfile.close()

mb2hmdb = {}

for line in mb2hmdblines:
	lsplit = line.split("\t")

	mb2hmdb[lsplit[0].strip().strip('"')] = lsplit[1].strip().strip('"')
	# print lsplit[0].strip().strip('"')

class EntryStruct:
	def __init__(self):
		self.accession = ''
		self.mol_name = ''
		self.formula = ''
		self.exact_mass = ''
		self.smiles = ''
		self.inchi = ''
		self.instrument = ''
		self.instrument_type = ''
		self.ms_type = ''
		self.ion_mode = ''
		self.collision_energy = ''
		self.ion_spray_voltage = ''
		self.ionization = ''
		self.retention_time = ''
		self.focused_ion_mz = ''
		self.focused_ion_type = ''
		self.data_procession = ''
		self.num_peaks = ''
		self.peaks = []
		self.normed_ints = []
		self.external_id = ''

	def printFields(self):
		print "ACCESSION:", self.accession
		print "NAME:", self.mol_name
		print "FORMULA:", self.formula
		print "EXACT MASS:", self.exact_mass
		print "SMILES:", self.smiles
		print "INCHI:", self.inchi
		print "INSTRUMENT:", self.instrument
		print "INSTRUMENT_TYPE:", self.instrument_type
		print "MS_TYPE:", self.ms_type
		print "ION_MODE:", self.ion_mode
		print "COLLISION_ENERGY:", self.collision_energy
		print "ION_SPRAY_VOLTAGE:", self.ion_spray_voltage
		print "IONIZATION:", self.ionization
		print "RETENTION TIME:", self.retention_time
		print "FOCUSED ION MZ:", self.focused_ion_mz
		print "FOCUSED ION TYPE:", self.focused_ion_type
		print "DATA PROCESSING:", self.data_procession
		print "NUMBER OF PEAKS:", self.num_peaks
		print "PEAKS:"
		# for i in self.peaks:
			# print i[0], i[1]
		print "// END //"

	def isValid(self):
		if (self.focused_ion_mz == ""):
			return False

		return True

mapNameToID = {}

dictMSTypes = {}

levels = {}

RT = 1.0

false_records = 0
no_prec_ion = 0

pos_records = 0
neg_records = 0

hmdb_found = 0
total_records = 0

for fname in filenames:
	mbfile = open(sys.argv[1] + fname.strip(), "r")
	mbfilelines = mbfile.readlines()
	mbfile.close()

	molname = ''
	entry = EntryStruct()

	peaks_open = False
	illegal_record = False
	first_name = True

	for mbline in mbfilelines:
		lsplit = mbline.strip().split()
		
		try:
			queryID = lsplit[0].strip()
		except:
			print "bla"

		if (queryID == "ACCESSION:"):
			entry.accession = mbline[11:].strip()

		if (queryID == "CH$NAME:"):
			if (first_name == True):
				entry.mol_name = mbline[9:].strip()
				first_name = False

		if (queryID == "CH$FORMULA:"):
			entry.formula = mbline[12:].strip()

		if (queryID == "CH$EXACT_MASS:"):
			entry.exact_mass = mbline[15:].strip()

		if (queryID == "CH$SMILES:"):
			entry.smiles = mbline[11:].strip()

		if (queryID == "CH$IUPAC:"):
			entry.inchi = mbline[10:].strip()

		if (queryID == "AC$INSTRUMENT:"):
			entry.instrument = mbline[14:].strip()

		if (queryID == "AC$INSTRUMENT_TYPE:"):
			entry.instrument_type = mbline[20:].strip()
			# found = typestring.find("LC-ESI")

			#if (found != -1):
			#	entry.meas_type = typestring

		if (queryID == "AC$MASS_SPECTROMETRY:"):
			if (lsplit[1].strip() == "MS_TYPE"):
				entry.ms_type = mbline[30:].strip()

			if (lsplit[1].strip() == "ION_MODE"):
				entry.ion_mode = mbline[31:].strip()

			if (lsplit[1].strip() == "COLLISION_ENERGY"):
				entry.collision_energy = mbline[39:].strip()

			if (lsplit[1].strip() == "ION_SPRAY_VOLTAGE"):
				entry.ion_spray_voltage = mbline[40:].strip()

			if (lsplit[1].strip() == "IONIZATION"):
				entry.ionization = mbline[33:].strip()

		if (queryID == "AC$CHROMATOGRAPHY:"):
			if (lsplit[1].strip() == "RETENTION_TIME"):
				entry.retention_time = mbline[34:].strip()

		if (queryID == "MS$FOCUSED_ION:"):
			if (lsplit[1].strip() == "PRECURSOR_M/Z"):
				entry.focused_ion_mz = mbline[30:].strip()

			if (lsplit[1].strip() == "FULL_SCAN_FRAGMENT_ION_PEAK"):
				entry.focused_ion_mz = mbline[44:].strip()

			if (lsplit[1].strip() == "PRECURSOR_TYPE"):
				entry.focused_ion_type = mbline[31:].strip()

			if (lsplit[1].strip() == "ION_TYPE"):
				entry.focused_ion_type = mbline[25:].strip()

		if (queryID == "CH$LINK:"):
			entry.external_id = lsplit[1].strip() + "_" + lsplit[2].strip()

		if (queryID == "PK$NUM_PEAK:"):
			entry.num_peaks = mbline[13:].strip()

		if (queryID == "PK$PEAK:"):
			if (not peaks_open):
				peaks_open = True
				continue

		if (peaks_open):
			d_split = mbline.split()

			if (mbline.strip() == '//'):
				peaks_open = False
				# do sanity check
				if (len(entry.peaks) != int(entry.num_peaks)):
					print "NUM_PEAKS and actual number of collected peaks differ!", len(entry.peaks), entry.num_peaks
					illegal_record=True
				continue

			if (len(d_split) == 3):
				# print "appending ", float(d_split[0]), float(d_split[1])
				# print "from ", d_split[0], d_split[1]
				# if (float(d_split[0]) > 5000.0):
				#	print "BIG: ", fname, d_split[0], d_split[1], d_split[2]
				entry.peaks.append((float(d_split[0].strip()), float(d_split[1].strip())))
			else:
				print "PEAK ENTRY FAULTY!"
				illegal_record = True
				false_records += 1


	if (not entry.isValid()):
		no_prec_ion += 1
		#entry.printFields()
		# print entry.inchi, "\t", entry.exact_mass,  "\t", entry.focused_ion_type, "\t", entry.ms_type, "\t", entry.collision_energy

		# try to repair broken entry
		add_ion = entry.focused_ion_type

		if (add_ion == "[M+H]+"):
			entry.focused_ion_mz = str(float(entry.exact_mass) + 1.007276)

		if (add_ion == "[M-H]-"):
			entry.focused_ion_mz = str(float(entry.exact_mass) - 1.007276)

		if (add_ion == "M+" or add_ion == "[M]+" or add_ion == "[M]+*"):
			if (("/q+1" in entry.inchi) or ("/p+1" in entry.inchi)):
				entry.focused_ion_mz = entry.exact_mass
			else:
				entry.focused_ion_mz = str(float(entry.exact_mass) + 1.007276)

	# entry.printFields()
	if (not illegal_record and entry.isValid()):
		# normalize intensities
		# get max
		int_max = -1
		for p in entry.peaks:
			if (p[1] > int_max):
				int_max = p[1]

		# print "max: ", int_max

		for p in entry.peaks:
			entry.normed_ints.append((p[1]/int_max)*100)

		# build spectra, one for each precursor
		mz_split = entry.focused_ion_mz.split("/")
		prectype_split = entry.focused_ion_type.split("/")

		if (len(mz_split) == 0):
			mz_split = entry.focused_ion_mz

		if (len(prectype_split) == 0):
			prectype_split = entry.focused_ion_type

		ch_i = 0
		charge_value = 1

		precursors = []

		for m in mz_split:
			prec = Precursor()
			prec.setMZ(float(m))

			# determine charge
			if (ch_i < len(prectype_split)):
				# retrieve charge from pseudomolecular ion
				ion_split = prectype_split[ch_i].split("]")

				if (len(ion_split) == 2):
					tmp_charge = ion_split[1].strip()

					charge_value = 1
					tmp_ch_split = tmp_charge.split("+")

					if (len(tmp_ch_split) == 2):
						if (tmp_ch_split[1] != "" and tmp_ch_split[1] != "*"):
							charge_value = int(tmp_ch_split[1])

					tmp_ch_split = tmp_charge.split("-")

					if (len(tmp_ch_split) == 2):
						if (tmp_ch_split[1] != "" and tmp_ch_split[1] != "*"):
							charge_value = int(tmp_ch_split[1])

						charge_value *= -1

					# print charge_value
					ch_i += 1

			prec.setCharge(charge_value)

			precursors.append(prec)

			spec = MSSpectrum()
			spec.setRT(RT)
			spec.setPrecursors(precursors)
			spec.setMetaValue("Massbank_Accession_ID", entry.accession)

			hmdb_id = 'NA'

			if (mb2hmdb.has_key(entry.accession)):
				hmdb_id = mb2hmdb[entry.accession]
				hmdb_found += 1


			spec.setMetaValue("HMDB_ID", hmdb_id)
			spec.setMetaValue("Inchi_String", entry.inchi)
			spec.setMetaValue("SMILES_String", entry.smiles)
			spec.setMetaValue("Precursor_Ion", entry.focused_ion_type)
			spec.setMetaValue("Sum_Formula", entry.formula)
			spec.setMetaValue("Precursor_Charge", charge_value)

			# check metabolite name for illegal < or > characters and replace them
			entry.mol_name = entry.mol_name.replace("<", "&lt;")
			# entry.mol_name = entry.mol_name.replace("<", "&lt;")


			spec.setMetaValue("Metabolite_Name", entry.mol_name)
			# print entry.mol_name, entry.focused_ion_mz

			spec.setType(SpectrumSettings().SpectrumType().PEAKS)

			type_split = entry.ms_type.split("S")

			ms_level = 1

			if (entry.collision_energy != ""):
				ms_level = 2

			if (len(type_split) > 1 and type_split[1] != ""):
				ms_level = int(type_split[1])
			#else:
			#	print entry.ms_type, ms_level, entry.collision_energy

			if (not levels.has_key(ms_level)):
				levels[ms_level] = 0

			levels[ms_level] += 1

			spec.setMSLevel(ms_level)

			# print ms_level, entry.focused_ion_mz, entry.focused_ion_type

			for p_i in range(0,len(entry.peaks)):
				peak = Peak1D()
				peak.setMZ(entry.peaks[p_i][0])
				peak.setIntensity(entry.normed_ints[p_i])
				spec.push_back(peak)

			msexp.addSpectrum(spec)
			RT += 1.0
			total_records += 1

mzmlfile.store("MBSpectra.mzML", msexp)

print false_records, no_prec_ion
print hmdb_found, total_records
