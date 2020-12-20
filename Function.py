from numbers import Number
from collections.abc import Callable
from functools import reduce
from operator import mul

import math

class Function:
    """An abstract class representing mathematical functions"""

    # Variable of the functions (used in string representation of the functions)
    # customisation might get added at some point
    variable = "x"

    def __init__(self, args=None):
        if args is not None:
            self.values = list(args)
        self.iszero = False

    def __call__(self, *args: Number) -> Number:
        return self.calculate_value(*args)

    def calculate_value(self, *args: Number) -> Number:
        # Calculates the functions value at specific point
        # (or points if multidimensional)
        raise NotImplementedError

    def __add__(self, other):
        if isinstance(other, Number):
            # adds support for summing numbers (constants) directly to functions 
            other = Constant(coef=other)
        elif isinstance(other, Sum):
            # triggers Sum's more clever addition
            return other + self
        return Sum(self, other)

    def __sub__(self, other):
        return self + (-1*other)

    def __mul__(self, other):
        if isinstance(other, Number):
            # support for multiplying by constants
            other = Constant(coef=other)
        elif isinstance(other, Product):
            # triggers Product multiplication
            return other * self
        return Product(self, other)

    def __matmul__(self, other):
        # The @ operator is used here to represent composition
        # i.e. outer @ inner
        return Composition(self, other)

    def __neg__(self):
        return self * (-1)

    def __truediv__(self, other):
        return self * Composition(Reciprocal(), other)

    def __radd__(self, other):
        return Constant(coef=other) + self

    def __rsub__(self, other):
        return Constant(coef=other) - self

    def __rmul__(self, other):
        return self * other

    def __rtruediv__(self, other):
        return other * Composition(Reciprocal(), self)

    def __pow__(self, other):
        if isinstance(other, Number):
            return power(degree=other) @ self
        return NotImplemented
    
    def derivative(self):
        raise NotImplementedError

    def integral(self):
        raise NotImplementedError


class SpecialFunction(Function):

    def __init__(self, kwargs={}):
        super().__init__()
        self.coef = kwargs.setdefault("coef", 1)
        self.kwargs = kwargs
        self.iszero = self.coef == 0

    def calculate_value(self, x):
        return self.func(x)

    # will currently cause errors because the functions take different input values
    # could be solved by functions taking parameters as **kwargs
    def __mul__(self, other):
        if isinstance(other, Number):
            kwargs = self.kwargs.copy()
            kwargs["coef"] *= other
            return type(self)(**kwargs)
        return super().__mul__(other)

    def coef_str(self):
        # returns the coef of the function as a string
        # + and - 1 are handled differently
        return str({1: "", -1: "-", 1.0: "", -1.0: "-"}.get(self.coef, str(self.coef) + "*"))

    def string(self, *args):
        # every function must 
        return NotImplemented

    def __str__(self):
        return self.string()


class Constant(SpecialFunction):
    """Returns a constant function"""
    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef
    
    def string(self, *args):
        return str(self.coef)

    def derivative(self):
        return Zero()

    def integral(self):
        return power(coef=self.coef)

class Zero(Constant):
    def __init__(self, **kwargs):
        kwargs["coef"] = 0
        super().__init__(**kwargs)

    def derivative(self):
        return self

def power(**kwargs):
    degree = kwargs.setdefault("degree", 1)
    if float(degree).is_integer():
        if degree > 0:
            return Monomial(**kwargs)
        elif degree < 0:
            return Reciprocal(**kwargs)
        else:
            return Constant(**kwargs)
    else:
        return Root(**kwargs)  

class Monomial(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.degree = kwargs.get("degree", 1)
        self.func = lambda x: self.coef * x**self.degree
    
    def string(self, var=Function.variable):
        # for example 3*x^2
        degree = "" if self.degree == 1 else f"^{self.degree}"
        return self.coef_str() + f"{var}{degree}"
    
    def derivative(self):
        if self.degree == 1:
            return Constant(coef=self.coef) 
        return self.degree * self.coef * power(degree=(self.degree - 1))

    def integral(self):
        degree = self.degree + 1
        return power(coef=self.coef/degree, degree=degree)

class Root(Monomial):

    def __init__(self, **kwargs):
        kwargs.setdefault("degree", 1/2)
        super().__init__(**kwargs)

    def string(self, var= "(" + Function.variable + ")"):
        # √x instead of 2√x
        n = "" if self.degree == 1/2 else 1/self.degree
        return self.coef_str() + f"{n}√{var}"

    def derivative(self):
        return self.coef * self.degree * Reciprocal(degree=1) @ Root(degree=pow(self.degree - 1, -1))


class Reciprocal(Monomial):

    def __init__(self, **kwargs):
        # without this calling Reciprocal() creates a 1st degree polynomial
        kwargs.setdefault("degree", -1)
        super().__init__(**kwargs)

    def string(self, var=Function.variable):
        degree = "" if self.degree == -1 else f"^{-self.degree}"
        return f"{self.coef}/{var}{degree}"
        
    def __mul__(self, other):
        if isinstance(other, Number):
            return Reciprocal(degree=-self.degree, coef=self.coef*other)
        else:
            return super().__mul__(other)

    def integral(self):
        if self.degree == -1:
            return Logarithm(coef=self.coef)
        return super().integral()

class Exponential(SpecialFunction):

    # adds support for "e" as a base
    to_numeric = {"e": math.e}

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.base = kwargs.get("base", "e")
        self.func = lambda x: self.coef*pow(Exponential.to_numeric.get(self.base), x)

    def string(self, var="(" + Function.variable + ")"):
        return self.coef_str() + f"{self.base}^{var}" 
    
    def derivative(self):
        return math.log(Exponential.to_numeric.get(self.base))*Exponential(base=self.base)

    # def __mul__(self, other):
    #     if isinstance(other, Number):
    #         return Exponential(base=self.base, coef=self.coef*other)
    #     else:
    #         return super().__mul__(other)

    def integral(self):
        kwargs = self.kwargs.copy()
        kwargs["coef"] /= math.log(Exponential.to_numeric[self.base])
        return Exponential(**kwargs)

class Logarithm(SpecialFunction):

    to_numeric = {"e": math.e}

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.base = kwargs.get("base", "e")
        self.func = lambda x: self.coef*math.log(x, Logarithm.to_numeric.get(self.base))

    def string(self, var="(" + Function.variable + ")") -> str:
        base = "ln" if self.base == "e" else f"log{self.base}"
        return self.coef_str() + f"{base}{var}"

    def derivative(self) -> Function:
        return self.coef/math.log(math.e)*Reciprocal(degree=1)

    def __mul__(self, other) -> Function:
        if isinstance(other, Number):
            return Logarithm(base=self.base, coef=self.coef*other)
        else:
            return super().__mul__(other)
    
class Sine(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef*math.sin(x)

    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str()  + f"sin{var}"

    def derivative(self) -> Function:
        return Cosine(coef=self.coef)

    def integral(self):
        return -Cosine(**self.kwargs)

class Cosine(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef*math.cos(x)
    
    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str()  + f"cos{var}"

    def derivative(self) -> Function:
        return -Sine(coef=self.coef)

    def integral(self):
        return Sine(**self.kwargs)

class Tangent(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef * math.tan(x)
    
    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str()  + f"tan{var}"

    def derivative(self) -> Function:
        return Tangent(coef=self.coef)**2 + 1

class Arcsin(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef * math.asin(x)

    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str() + f"arcsin{var}"

    def derivative(self) -> Function:
        return self.coef / (power(degree=1/2) @ (1 - power(degree=2)))

class Arccos(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef * math.acos(x)

    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str() + f"arccos{var}"

    def derivative(self) -> Function:
        return -self.coef / (power(degree=1/2) @ (1 - power(degree=2)))

class Arctan(SpecialFunction):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        self.func = lambda x: self.coef * math.atan(x)

    def string(self, var="(" + Function.variable + ")") -> str:
        return self.coef_str() + f"arctan{var}"

    def derivative(self) -> Function:
        return self.coef / (1 + power(degree=2))


class Sum(Function):

    def __init__(self, *args):
        super().__init__(args)
        # filters out values that are zero
        # note that self.values has already been instantiniated here
        self.values = list(filter(lambda f: not f.iszero, self.values))
        # zero if sum is empty
        self.iszero = self.values == []

    def calculate_value(self, *args: Number) -> Number:
        return sum(f(args[0]) for f in self.values)

    def __str__(self) -> str:
        # returns 'f(x) + g(x) + h(x)'
        return " + ".join(map(lambda f: f.string(), self.values))

    def string(self) -> str:
        # this will be called when the sum is inside another
        # function, like a product
        return "(" + str(self) + ")"

    def __add__(self, other: Function) -> Function:
        # without this a + b + c would become 
        # Sum(Sum(a, b), c) instead of Sum(a, b, c)
        if isinstance(other, Sum):
            return Sum(*(self.values + other.values))
        return Sum(*self.values, other)

    def __mul__(self, other: Function) -> Function:
        if isinstance(other, Number):
            return Sum(map(lambda f: other*f, self.values))
        else:
            return super().__mul__(other)

    def derivative(self) -> Function:
        # all zero valued functions can be removed from the sum
        nonzero = filter(lambda f: not f.iszero, self.values)
        return Sum(*map(lambda f: f.derivative(), nonzero))

    def simplify(self) -> Function:
        S = Sum()
        funcs = {}
        # sorts the functions to compatible sets
        for func in self.values:
            # every function has a (non-unique) id, which tells which functions
            # it is compatible with
            # for example x^2 and 2*x^2 have the same, so they can be summed
            comp_id = func.compatibility_id
            if comp_id not in funcs.keys():
                funcs[comp_id] = []
            funcs[comp_id] += func
        
        for L in funcs.items():
            if len(L) >= 2:
                S += reduce(lambda x, y: x.comp_sum(y), L[1:], L[0])
            else: S += L[0]
        
        return S


    def integral(self):
        return Sum(*map(lambda f: f.integral(), self.values))


class Product(Function):

    def __init__(self, *args):
        super().__init__(args)
        # 0 if any factors are 0
        self.iszero =  any(f.iszero for f in self.values)

    @staticmethod
    def product(init: Number, *args: Number) -> Number:
        # a help function for calculating the products 
        return reduce(mul, *args, init)

    def calculate_value(self, *args: Number) -> Number:
        return Product.product(map(lambda f: f(args[0]), self.values))

    def __str__(self) -> str:
        # i.e. a * b * c
        return " * ".join(map(lambda f: f.string(), self.values))

    def string(self) -> str:
        return str(self)

    def __mul__(self, other: Function) -> Function:
        # because of this a * b * c is 
        # Product(a, b, c) not Product(Product(a, b), c)
        if isinstance(other, Product):
            return Product(*(self.values + other.values))
        return Product(*self.values, other)

    def derivative(self) -> Function:
        if self.iszero:
            return Zero()
        return Sum(*map(lambda val: Product(*filter(lambda v: v != val, self.values), val.derivative()), self.values))


class Composition(Function):

    def __init__(self, outer, inner):
        self.values = tuple((outer, inner))
        self.__outer = outer
        self.__inner = Zero() if inner.iszero else inner
        # the function is zero if the outer function is zero
        self.iszero = self.__outer.iszero

    def calculate_value(self, *args: Number) -> Number:
        return self.__outer(self.__inner(args[0]))

    def __str__(self) -> str:
        s = "(" + str(self.__inner) + ")"
        return f"{self.__outer.string(s)}"

    def string(self, var="") -> str:
        return str(self)

    def derivative(self) -> Function:
        return self.__inner.derivative() * (self.__outer.derivative() @ self.__inner)


f = Cosine() + Exponential()
print(f)
print(f.derivative())
print(f.integral())