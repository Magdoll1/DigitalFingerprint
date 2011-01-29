import os
import sys
from DF import DF, DFReader, DFWriter
"""
Quick script to combine different date runs of the same person
ex: V1.CD-8_090127 and V1.CD-8_090109 into the same sample
"""

#df_filename = 'Qin_MetaHit.justSpanish_dateseparate.DF'
df_filename = 'Qin_MetaHit.justDanish_dateseparate.DF'

df_by_person = {}
for df in DFReader(open(df_filename)):
	person = df.name[:df.name.find('_')]
	if person in df_by_person:
		# combine the two!
		print >> sys.stderr, "combining {0} with {1}".format(df_by_person[person].name, df.name)
		df_by_person[person] += df
		df_by_person[person].name = person
	else:
		df_by_person[person] = df


#f = open('Qin_MetaHit.justSpanish_datecombined.DF', 'w')
f = open('Qin_MetaHit.justDanish_datecombined.DF', 'w')
w = DFWriter(f)
w.writes(df_by_person.itervalues())
f.close()
