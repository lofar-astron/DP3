import os
import sys

name, _ = os.path.splitext(sys.argv[1])
with open(sys.argv[1]) as f:
	lines = f.readlines()

with open('schemas/' + sys.argv[2] + '.yml', 'w') as w:
	w.write('description: >-\n')
	w.write('  ""\n')
	w.write('inputs:\n')
	for line in lines:
		print(line)
		[_, parameter, the_type, the_default, doc, _] = line.replace("''|''", "''OROR''").replace("''||''", "''OROROR''").split('|')

		w.write('  ' + parameter.replace('<step>.', '').strip().replace('.', '&#46;') + ':\n')

		if len(the_default.strip()) > 0:
			w.write('    default: ' + the_default.replace("''", '"\\"\\""').replace('""', '"\\"\\""').replace('[]', '"[]"').strip() + '\n')

		w.write('    type: ' + the_type.replace(' vector', '?').strip() + '\n')

		doc = doc.strip()
		if doc.endswith('.'):
			doc = doc[:-1] + ' `.`'
		elif doc.endswith('?'):
			doc = doc[:-1] + ' `?`'
		else:
			doc = doc + ' `.`'
		w.write('    doc: >-\n')
		w.write('      ' + doc.replace('OROR', '|').replace('OROROR', '||').replace('\\\\ ', '').replace("''", '``') + '\n')
