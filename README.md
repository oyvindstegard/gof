# Game of Life

This is an HTML5 canvas experiment implementing J.H. Conway's famous "Game of
Life" cellular automaton game. In its pure mathematical form, the game requires
an infinite grid. But given the constraint of finite computing resources, this
would be impractical, so the universe is limited in size. This may affect
evolutions to not be 100% accurate mathematically, depending on the amount of
grid space required. It's still fun, though.

## About the code

The engine code `gof.js` is plain JavaScript from 2014. The canvas drawing
routines and engine could probably be heavily optimized or re-written for better
performance, but I currently do not have the time to work on it. The user
interface `index.html` requires jQuery.

## See it running

Here: https://stegard.net/game-of-life/


