import sys

ZhuYintoBig5 = {}

with open(sys.argv[1], 'r', encoding = 'big5-hkscs') as f:
	for line in f:
		[word, zhuyin] = line.strip('\n').split(' ')
		
		if '/' in zhuyin:
			zys = zhuyin.strip().split('/')
			for zy in zys:
				if zy[0] not in ZhuYintoBig5:
					ZhuYintoBig5[zy[0]] = [word]
				else:
					if word not in ZhuYintoBig5[zy[0]]:
						ZhuYintoBig5[zy[0]].append(word)
					else:
						None
		else:
			zhuyin = zhuyin.strip()
			if zhuyin[0] not in ZhuYintoBig5:
				ZhuYintoBig5[zhuyin[0]] = [word]
			else:
				if word not in ZhuYintoBig5[zhuyin[0]]:
					ZhuYintoBig5[zhuyin[0]].append(word)
				else:
					None

with open(sys.argv[2], 'w', encoding = 'big5-hkscs') as f:
	for (zhuyin, words) in ZhuYintoBig5.items():
		f.write(zhuyin+'	'+' '.join(words)+'\n')
		print(zhuyin+'	'+' '.join(words)+'\n')
		for word in words:
			f.write(word+'	'+word+'\n')
			print(word+'	'+word+'\n')
					


