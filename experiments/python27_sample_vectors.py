"""Run under actual CPython2.7 only; explicit frozen CDF and source draw order."""
import json,random,sys
j=json.load(open(sys.argv[1]));out=[]
for seed in j['seeds']:
 random.seed(seed);u=random.random();degree=next(i+1 for i,c in enumerate(j['cdf']) if c>u)
 out.append({'seed':seed,'first_uniform':u,'degree':degree,'indices':random.sample(xrange(j['K']),degree)})
json.dump({'python_version':sys.version,'vectors':out},open(sys.argv[2],'w'),sort_keys=True,separators=(',',':'))
