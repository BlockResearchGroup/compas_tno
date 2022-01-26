import matlab.engine

# future = matlab.engine.start_matlab(background=True)
future = matlab.engine.connect_matlab(background=True)
eng = future.result()
eng.cvx_begin(nargout=0)
print(matlab.engine.find_matlab())
