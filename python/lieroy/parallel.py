
import importlib

class FunctionWrapper():
    def __init__(self, function_name, module_name):
        self.function_name = function_name
        self.module_name = module_name

    def __call__(self, *args, **kwargs):
        module = importlib.import_module(self.module_name)
        return module.__dict__[self.function_name](*args, *kwargs)
