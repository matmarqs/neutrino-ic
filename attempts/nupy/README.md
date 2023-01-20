# nupy
nusim but in python

## atualizações

```sh
time python3 ivp.py > output.txt
```

A probabilidade de sobrevivência está sempre travada no 1, para pelo menos
os raios menores que 0.942. A partir daí o código para de funcionar.

Ele emite o seguinte erro:

```sh
ValueError: A value in x_new is above the interpolation range.
```
