import matlab.engine

# Run the following lines in the terminal if you have trouble

future = matlab.engine.connect_matlab(background=True)
eng = future.result()
eng.cvx_begin(nargout=0)
print(matlab.engine.find_matlab())
