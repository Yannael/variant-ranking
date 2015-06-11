data = xrange(1, 10001)
xrangeRDD = sc.parallelize(data, 8)


(sc
 .parallelize(data)
 .map(lambda y: y - 1)
 .filter(lambda x: x < 10)
 .collect())
 
 
 wordsList = ['cat', 'elephant', 'rat', 'rat', 'cat']
wordsRDD = sc.parallelize(wordsList, 4)
