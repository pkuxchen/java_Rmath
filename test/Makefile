JAVA_HOME = ..

CLASS = rmathtest.class
RUN = rmath_test_result.txt

all: $(CLASS) $(RUN)

$(CLASS):
	javac -classpath $(JAVA_HOME)/class:. rmathtest.java 

$(RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:. rmathtest > rmath_test_result.txt

clean: 
	rm -f *.txt *.class
