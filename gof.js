// Game of Life in JavaScript - version 1-beta.
// Copyright (C) 2014 Øyvind Stegard.

'use strict';

// TODO touch event handling (for drawing on and resizing grid)

// GOF namespace object
var GOF = {};

// Main game code
GOF.Game = function(canvas, options) {
    var that = this; // pointer to object instance for private methods

    if (typeof canvas == "string") {
        canvas = document.getElementById(canvas);
    }
    
    if (! (canvas && canvas.getContext)) {
        throw new Error("Canvas not supported or element not found");
    }

    // Grid dimensions, cell model structures and drawing objects
    var cols, rows, visibleCols, visibleRows;
    var visibleOffsetCols;
    var visibleOffsetRows;
    var cells;
    var iteration;
    var nextIterationStartRow = -1;
    var nextIterationEndRow = -1;
    var ctx;
    var shapeCache = null; // internal canvas for caching cell shapes as bitmaps
    var initialize = function(initialPattern) {
        visibleCols = Math.floor((canvas.width - that.padding)
                                          / (that.padding + that.cellSize));
        visibleRows = Math.floor((canvas.height - that.padding)
                                          / (that.padding + that.cellSize));

        // Prefer even numbered visible cols and rows
        visibleCols -= visibleCols % 2;
        visibleRows -= visibleRows % 2;

        initialPattern = initialPattern ? initialPattern : that.initialPattern;

        var minSizeCols = initialPattern.sizeCols ? initialPattern.sizeCols : 3;
        var minSizeRows = initialPattern.sizeRows ? initialPattern.sizeRows : 3;
        var gridSize = calculateGridSize(visibleCols, visibleRows, minSizeCols, minSizeRows);
        var currentCols = cols;
        var currentRows = rows;
        cols = gridSize.cols;
        rows = gridSize.rows;
        visibleOffsetCols = gridSize.visibleOffsetCols;
        visibleOffsetRows = gridSize.visibleOffsetRows;

        iteration = 0;
        nextIterationStartRow = -1;
        nextIterationEndRow = -1;

        if (cells) {
            cells = resizeCellGrid(cells, currentCols, currentRows, cols, rows);
        } else {
            cells = createCellGrid(cols, rows, visibleCols, visibleRows, initialPattern);
        }

        shapeCache = null;
        ctx = canvas.getContext("2d");
        updateCanvas();
    };

    // Returns total grid size and visibility offsets as an object with
    // properties "cols", "rows", "visibleOffsetCols" and "visibleOffsetRows".
    var calculateGridSize = function(visibleCols, visibleRows, minSizeCols, minSizeRows) {
        var cols = Math.max(visibleCols, minSizeCols);
        cols += Math.max(100-cols, 4);
        var rows = Math.max(visibleRows, minSizeRows);
        rows += Math.max(100-rows, 4);

        var visibleOffsetCols = Math.floor((cols-visibleCols)/2);
        var visibleOffsetRows = Math.floor((rows-visibleRows)/2);
        return { "cols": cols,
                 "rows": rows,
                 "visibleOffsetCols": visibleOffsetCols,
                 "visibleOffsetRows": visibleOffsetRows };
    };
    
    /* Initialize a new cell array with optional predicate
       function for which cells should be set to "alive" initially
    */
    var createCellGrid = function(cols, rows, visibleCols, visibleRows, initialPattern) {
        var cells = new Array(cols * rows);
        var pfunc = initialPattern && initialPattern["pfunc"] ? initialPattern["pfunc"] : null;
        for (var i=0; i<cells.length; i++) {
            cells[i] = { "state": false, "otherState":false};
            if (pfunc) {
                cells[i].state = pfunc(i % cols, Math.floor(i / cols), cols, rows, visibleCols, visibleRows);
            }
        }
        return cells;
    };

    /* Return a resized cell grid based on cell states in old grid.
       Anchor point for grid resize is at the center. Keep cols/rows to
       even numbers to avoid center drift when sizing up/down.
    */
    var resizeCellGrid = function(cells, cols, rows, newCols, newRows) {
        if (cols == newCols && rows == newRows) {
            return cells;
        }
        
        if (cols <= newCols) {
            var copyColStart = 0;
            var copyColDest = Math.floor(newCols/2 - cols/2);
            var copyColDestEnd = copyColDest + cols;
        } else {
            var copyColStart = Math.floor((cols - newCols)/2);
            var copyColDest = 0;
            var copyColDestEnd = copyColDest + newCols;
        }
        if (rows <= newRows) {
            var copyRowStart = 0;
            var copyRowDest = Math.floor(newRows/2 - rows/2);
            var copyRowDestEnd = copyRowDest + rows;
        } else {
            var copyRowStart = Math.floor((rows - newRows)/2);
            var copyRowDest = 0;
            var copyRowDestEnd = copyRowDest + newRows;
        }

        var newSize = newCols * newRows;
        var newCells = new Array(newSize);
        var row = -1;
        for (var i=0; i<newSize; i++) {
            var col = i % newCols;
            if (col == 0) ++row;
            if (col >= copyColDest && col < copyColDestEnd
                 && row >= copyRowDest && row < copyRowDestEnd) {
                var sourceCol = copyColStart + col - copyColDest;
                var sourceRow = copyRowStart + row - copyRowDest;
                newCells[i] = cells[sourceRow*cols + sourceCol];
            } else {
                newCells[i] = {"state":false, "otherState":false};
            }
        }
        
        return newCells;
    };
    
    /* Determine next state for a single cell, based on current state of cell
       and number of live neighbours. Returns true if cell at index idx should be
       alive for the next iteration, false otherwise.

       The rules for Game of Life are:
         1. Any live cell with fewer than two live neighbours dies, as if caused by under-population.
         2. Any live cell with two or three live neighbours lives on to the next generation.
         3. Any live cell with more than three live neighbours dies, as if by overcrowding.
         4. Any dead cell with exactly three live neighbours becomes a live cell, as if by reproduction.

       Source: http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life#Rules
    */
    var nextCellState = function(cells, i) {
        var nCount = 0;
        for (var rowbase=(i-cols < 0 ? i : i-cols); rowbase <= i+cols && rowbase < cells.length; rowbase+=cols) {
            if (rowbase != i) {
                if (cells[rowbase].state) ++nCount;
            }
            var col = rowbase % cols;
            if (col == 0) {
                if (cells[rowbase+1].state) ++nCount;
            } else if (col == cols - 1) {
                if (cells[rowbase-1].state) ++nCount;
            } else {
                if (cells[rowbase+1].state) ++nCount;
                if (cells[rowbase-1].state) ++nCount;
            }
        }
        
        return cells[i].state ? (nCount == 2 || nCount == 3) : (nCount == 3);
    };

    /* Do canvas update according to cell grid state, optionally incremental redraw.
       If 'incremental' is true, then all cells whose state is different from previous
       state are redrawn. If it's a number, then it's interpreted as index of a single
       cell to redraw.
    */
    var updateCanvas = function(incremental) {
        var cellIndex = null;
        if (typeof incremental == "number") {
            cellIndex = incremental;
            incremental = true;
        }
        
        if (!incremental) {
            // Draw background
            ctx.fillStyle = getBackgroundFillStyle();
            ctx.fillRect (0, 0, canvas.width, canvas.height);
        }

        var shape = that.theme.cellShape == "circle" ? "circle" : "square";
        if (!incremental && shape != "circle" && that.padding == 0) {
            ctx.fillStyle = that.theme.cellColorDead;
            var gridCoords = gridCanvasCoords();
            ctx.fillRect(gridCoords[0], gridCoords[1], gridCoords[2], gridCoords[3]);
        }
        
        // Live cells
        var start = cellIndex ? cellIndex : 0;
        var end = cellIndex ? cellIndex + 1 : cells.length;
        ctx.fillStyle = that.theme.cellColorAlive;
        for (var i=start; i<end; i++) {
            if (cells[i].state) {
                if (incremental && cells[i].otherState) continue;
                var cellCoords = calculateCellCoords(i);
                if (cellCoords) {
                    fillCellShape(ctx, cellCoords, shape, true);
                }
            }
        }
        
        // Dead cells
        ctx.fillStyle = that.theme.cellColorDead;
        if (incremental || that.padding > 0 || shape == "circle") {
            for (var i=start; i<end; i++) {
                if (!cells[i].state) {
                    if (incremental && !cells[i].otherState) continue;
                    var cellCoords = calculateCellCoords(i);
                    if (cellCoords) {
                        fillCellShape(ctx, cellCoords, shape, false);
                    }
                }
            }
        }
    };

    // Fill cell shape. Fill style should be set prior to calling this function.
    // Caches rendered circle shapes for faster drawing.
    var fillCellShape = function(ctx, cellCoords, shape, state) {
        var cellX = cellCoords[0];
        var cellY = cellCoords[1];
        var cellWidth = cellCoords[2];
        var cellHeight = cellCoords[3];
        if (shape == "circle") {
            // Need background-clear for circle redraws
            ctx.save();
            ctx.fillStyle = getBackgroundFillStyle();
            ctx.fillRect(cellX, cellY, cellWidth, cellHeight);
            ctx.restore();
            if (!shapeCache) {
                shapeCache = document.createElement("canvas");
                shapeCache.width = cellWidth * 2;
                shapeCache.height = cellHeight;
                var cacheCtx = shapeCache.getContext("2d");
                var x=0, y=0;
                var radius = Math.floor(cellWidth/2);
                var centerX = radius + x;
                var centerY = radius + y;
                cacheCtx.fillStyle = that.theme.cellColorAlive;
                cacheCtx.beginPath();
                cacheCtx.moveTo(x, y);
                cacheCtx.arc(centerX, centerY, radius, 0, 2*Math.PI, false);
                cacheCtx.closePath();
                cacheCtx.fill();
                x += cellWidth;
                centerX = radius + x;
                cacheCtx.fillStyle = that.theme.cellColorDead;
                cacheCtx.beginPath();
                cacheCtx.moveTo(x, y);
                cacheCtx.arc(centerX, centerY, radius, 0, 2*Math.PI, false);
                cacheCtx.closePath();
                cacheCtx.fill();
            }
           
            if (state) {
                ctx.drawImage(shapeCache, 0, 0, cellWidth, cellHeight, cellX, cellY, cellWidth, cellHeight);
            } else {
                ctx.drawImage(shapeCache, cellWidth, 0, cellWidth, cellHeight, cellX, cellY, cellWidth, cellHeight);
            }

        } else { // if (shape == "circle") { ...
            ctx.fillRect(cellX, cellY, cellWidth, cellHeight);
        }
    };

    // Returns fill style for background based on theme
    var getBackgroundFillStyle = function() {
        if (that.theme.bgGradient) {
            return that.theme.bgGradient;
        }
        if (that.theme.bgColorEnd) {
            var gradient = ctx.createLinearGradient(0,0,0,canvas.height);
            gradient.addColorStop(0, that.theme.bgColor);
            gradient.addColorStop(1, that.theme.bgColorEnd);
            that.theme.bgGradient = gradient;
            return gradient;
        } else {
            return that.theme.bgColor;
        }
    }

    // Returns an array with the first two elements being the top-left coordinates
    // where grid should be drawn on canvas, and the 3rd and 4th elements being
    // the grid pixel size in x and y direction. 
    var gridCanvasCoords = function() {
        var sizeX = (that.cellSize + that.padding) * visibleCols - that.padding;
        var originX = Math.floor((canvas.width - sizeX) / 2);
        var sizeY = (that.cellSize + that.padding) * visibleRows - that.padding;
        var originY = Math.floor((canvas.height - sizeY) / 2);
        return [originX, originY, sizeX, sizeY];
    };

    // Calculate canvas x,y (top-left) coordinates for cell at index i
    // Returns an array of length 4 with x coordinate at index 0, y coordinate at index 1,
    // size in x direction at index 2 and size in y direction at index 3
    var calculateCellCoords = function(i) {
        var col = i % cols;
        var row = Math.floor(i / cols);
        // Only draw visible area of grid
        if (col < visibleOffsetCols || col >= visibleCols + visibleOffsetCols) {
            return undefined;
        }
        if (row < visibleOffsetRows || row >= visibleRows + visibleOffsetRows) {
            return undefined;
        }
        col -= visibleOffsetCols;
        row -= visibleOffsetRows;
        var gridCoords = gridCanvasCoords();
        return [gridCoords[0] + (col % visibleCols)*(that.cellSize + that.padding),
                gridCoords[1] + row*(that.cellSize + that.padding),
                that.cellSize, that.cellSize];
    };

    // Calculate cell index from a canvas point
    // Returns cell index if point is on a visible grid cell, or undefined otherwise
    // (padding or outside visible viewport)
    var cellIndexFromCoords = function(x, y) {
        var gridCoords = gridCanvasCoords();
        x -= gridCoords[0];
        y -= gridCoords[1];
        var inCellX = x % (that.cellSize + that.padding);
        if (inCellX >= that.cellSize + (that.padding <= 5 ? that.padding : 0)) {
            return undefined;
        }
        var inCellY = y % (that.cellSize + that.padding);
        if (inCellY >= that.cellSize + (that.padding <= 5 ? that.padding : 0)) {
            return undefined;
        }
        var visibleCol = Math.floor(x / (that.cellSize + that.padding));
        var visibleRow = Math.floor(y / (that.cellSize + that.padding));

        if (visibleCol >= visibleCols || visibleRow >= visibleRows) {
            return undefined;
        }

        return (visibleRow+visibleOffsetRows)*cols + visibleCol + visibleOffsetCols;
    };

    // Get top-left mouse position for mouse event relative to canvas bounding box
    var mousePositionCanvas = function(mouseEvent) {
        var rect = canvas.getBoundingClientRect();
        return [mouseEvent.clientX - rect.left, mouseEvent.clientY - rect.top];
    };

    // Handle certain mouse events on canvas element
    var mouseEventHandler = function(event) {
        var button = 0;
        if ("buttons" in event) {
            button = event.buttons;
        } else if ("which" in event) {
            button = event.which;
        } else {
            button = event.button;
        }
        if (button) {
            var create = (button === 1);
            var pos = mousePositionCanvas(event);
            var i = cellIndexFromCoords(pos[0], pos[1]);
            if (i) {
                cells[i].otherState = cells[i].state;
                cells[i].state = create ? true : false;
                var row = Math.floor(i/cols);
                nextIterationStartRow = Math.min(nextIterationStartRow, row-1);
                nextIterationEndRow = Math.max(nextIterationEndRow, row+1);
                updateCanvas(i);
            }
        }
    };

    /* Public members */

    /*  Proceed to next iteration.
     *  Returns 'true' if there were cell changes, 'false' otherwise. This may
     *  be used to halt game if no more changes are occuring.
     */
    this.nextIteration = function() {
        var len = cells.length;
        var start = Math.max(0, nextIterationStartRow*cols);
        var end = nextIterationEndRow > 0 ? Math.min(nextIterationEndRow*cols + cols, len) : len;

        for (var i=start; i<end; i++) {
            cells[i].otherState = nextCellState(cells, i);
        }

        var changes = false;
        nextIterationStartRow = -1;
        nextIterationEndRow = -1;
        for (var i=start; i<end; i++) {
            var cell = cells[i];
            if (cell.state != cell.otherState) {
                changes = true;
                var s = cell.state;
                cell.state = cell.otherState;
                cell.otherState = s;
            }

            // Simple optimization to avoid computing cell states for dead areas at the top
            // and bottom of the grid.
            if (cell.state) {
                var row = Math.floor(i / cols);
                if (nextIterationStartRow == -1) {
                    nextIterationStartRow = row - 1;
                }
                nextIterationEndRow = row + 1;
            }
        }

        updateCanvas(true);
        ++iteration;
        return changes;
    };

    /* Various options */
    if (options && typeof options["theme"] == "object") {
        this.theme = options["theme"];
        // Inherit from default theme if missing essential bits
        if (! this.theme.bgColor) {
            this.theme.bgColor = GOF.themes["Default"].bgColor;
        }
        if (! this.theme.cellColorAlive) {
            this.theme.cellColorAlive = GOF.themes["Default"].cellColorAlive;
        }
        if (! this.theme.cellColorDead) {
            this.theme.cellColorDead = GOF.themes["Default"].cellColorDead;
        }
    } else {
        this.theme = GOF.themes["Default"];
    }
    this.cellSize = 8;
    this.padding = 2;
    if (typeof this.theme.cellSize == "number") {
        this.cellSize = this.theme.cellSize;
    }
    if (typeof this.theme.padding == "number") {
        this.padding = this.theme.padding;
    }
    
    // Specific cell size and padding can override theme defaults
    if (options && typeof options["cellSize"] == "number") {
        this.cellSize = options["cellSize"];
    } 
    if (options && typeof options["padding"] == "number") {
        this.padding = options["padding"];
    }
    this.initialPattern = GOF.initialPatterns["Default"];
    if (options && typeof options["initialPattern"] == "object") {
        this.initialPattern = options["initialPattern"];
    }
    /* Set size of one cell, and optionally the padding between each cell.
       Values are in canvas pixel units. Grid size will be set accordingly.
    */
    this.setCellSize = function(size, padding) {
        if (typeof size !== "number") {
            throw new Error("Cell size must be a number");
        }
        this.cellSize = Math.max(Math.floor(size),1);
        if (typeof padding == "number") {
            this.padding = Math.max(Math.floor(padding),0);
        }
        initialize();
    };
    this.clear = function() {
        cells = null;
        initialize(function() { return false;});
    }
    this.reset = function() {
        cells = null;
        initialize();
    }
    this.setTheme = function(theme) {
        if (typeof theme == "string") {
            theme = GOF.themes[theme];
        }
        if (theme) {
            this.theme = theme;
            initialize();
        }
    };
    this.setInitialPattern = function(initP) {
        this.initialPattern = initP;
    };

    // Used for getting array literal representation of current visible part of cell grid.
    this.getVisibleCellsAsArrayString = function() {
        if (!cells) return undefined;
        var a = "cols=" + cols + ", rows=" + rows + ", vcols=" + visibleCols + ", vRows=" + visibleRows + ", cells=\n[";
        for (var i=0; i<cells.length; i++) {
            var col = i % cols;
            var row = Math.floor(i/cols);
            if (col >= visibleOffsetCols && col < cols - visibleOffsetCols) {
                if (row >= visibleOffsetRows && row < rows - visibleOffsetRows) {
                    a += (cells[i].state ? "1" : "0");
                    if ((col+1) < cols - visibleOffsetCols || row < rows - visibleOffsetRows - 1) {
                        a += ","
                    }
                    if (col == cols - visibleOffsetCols - 1 && row < rows - visibleOffsetRows - 1) {
                        a += "\n";
                    }
                }
            }
        }
        return a + "]";
    };

    // Create initial cell grid, draw and set up canvas event handler
    initialize();
    canvas.addEventListener("mousemove", mouseEventHandler);
    canvas.addEventListener("mousedown", mouseEventHandler);

}; // ctor Game

/* Get position object for placement in grid. Positioning is kept within visible
   part of grid if there is enough space. Cols/rows may be outside grid if rectangle
   size is greater than grid size.
*/ 
GOF.position = function(cols, rows, vCols, vRows, sizeCols, sizeRows, pSpec, offsetCol, offsetRow) {
    if (typeof pSpec != "string") {
        pSpec = "center";
    }
    var colPosition;
    var rowPosition;
    if (pSpec.indexOf("-") == -1) {
        if (pSpec == "top" || pSpec == "bottom") {
            colPosition="center";
            rowPosition=pSpec;
        } else if (pSpec == "left" || pSpec == "right") {
            colPosition=pSpec;
            rowPosition="center";
        } else {
            colPosition="center";
            rowPosition="center";
        }
    } else {
        colPosition=pSpec.substring(0,pSpec.indexOf("-"));
        rowPosition=pSpec.substring(pSpec.indexOf("-")+1, pSpec.length);
    }

    var p = { "tlCol" : null, "brCol" : null, "tlRow" : null, "brRow" : null,
              "cols" : cols, "rows" : rows, "vCols" : vCols, "vRows" : vRows,
              "sizeCols" : sizeCols, "sizeRows" : sizeRows };
    if (colPosition == "left") {
        p.tlCol = Math.floor((cols-vCols)/2);
        p.brCol = p.tlCol + sizeCols;
    } else if (colPosition == "right") {
        p.tlCol = Math.floor((cols-vCols)/2) + vCols - sizeCols;
        p.brCol = p.tlCol + sizeCols;
    } else {
        p.tlCol = Math.floor(cols/2 - sizeCols/2);
        p.brCol = p.tlCol + sizeCols;
    }
    if (rowPosition == "top") {
        p.tlRow = Math.floor((rows-vRows)/2);
        p.brRow = p.tlRow + sizeRows;
    } else if (rowPosition == "bottom") {
        p.tlRow = Math.floor((rows-vRows)/2) + vRows - sizeRows;
        p.brRow = p.tlRow + sizeRows;
    } else {
        p.tlRow = Math.floor(rows/2 - sizeRows/2);
        p.brRow = p.tlRow + sizeRows;
    }
    if (typeof offsetCol == "number") {
        p.tlCol += Math.floor(offsetCol);
        p.brCol += Math.floor(offsetCol);
    }
    if (typeof offsetRow == "number") {
        p.tlRow += Math.floor(offsetRow);
        p.brRow += Math.floor(offsetRow);
    }

    return p;
};

/* Returns an array of cells representing a circle of a given radius (in cells).
   The returned array will have cols equal to (radius*2 + 1), and the same number for rows. */
GOF.circleCellsPredicate = function(radius, pSpec, offsetCol, offsetRow) {
    var pos = null;
    var cells = null;
    return function(i,j,cols,rows,vCols,vRows) {
        if (!pos || !(pos.cols == cols && pos.rows == rows && pos.vCols == vCols && pos.vRows == vRows)) {
            // Recompute
            var circleCols = radius*2 + 1;
            pos = GOF.position(cols, rows, vCols, vRows, circleCols, circleCols, pSpec, offsetCol, offsetRow);
            cells = new Array(circleCols*circleCols);
            for (var i=0; i<cells.length; i++) {
                cells[i] = 0;
            }
            var circleCell = function(x,y) {
                var shift = radius;
                cells[x+shift  + (y+shift)*circleCols] = 1;
                // Circle symmetry points follow:
                cells[x+shift  + (-y+shift)*circleCols] = 1;
                cells[-x+shift + (y+shift)*circleCols] = 1;
                cells[-x+shift + (-y+shift)*circleCols] = 1;
                cells[y+shift  + (x+shift)*circleCols] = 1;
                cells[-y+shift + (x+shift)*circleCols] = 1;
                cells[y+shift +  (-x+shift)*circleCols] = 1;
                cells[-y+shift + (-x+shift)*circleCols] = 1;
            };

            // Calculate live cells using Bresenham circle rasterizing algorithm
            var x = 0;
            var y = radius;
            var dE = 3 - 2*radius; // Initial error from circle line
            do {
                if (dE < 0) {
                    // Move to the right (default), update error
                    dE = dE + 4*x + 6;
                } else {
                    // Move diagonally and update error
                    dE = dE + 4*x - 4*(y--) + 10;
                }
                circleCell(x++, y);
            } while (x <= y); // Second octant of circle
        }
        
        if (i >= pos.tlCol && i < pos.brCol && j >= pos.tlRow && j < pos.brRow) {
            var cellIndex = (j-pos.tlRow)*pos.sizeCols + i - pos.tlCol;
            return (cells[cellIndex] == 1);
        }

        return false;
    }
};

/* Creates a predicate function based on array of cells. The 'cells' will
   be placed in grid according to positioning parameters.
 */
GOF.cellArrayPredicate = function(cells, sizeCols, sizeRows, pSpec, offsetCol, offsetRow) {
    var pos = null;
    return function(i, j, cols, rows, vCols, vRows) {
        if (!pos || !(pos.cols == cols && pos.rows == rows && pos.vCols == vCols && pos.vRows == vRows)) {
            // Recompute positioning
            pos = GOF.position(cols, rows, vCols, vRows, sizeCols, sizeRows, pSpec, offsetCol, offsetRow);
        }
        
        if (i >= pos.tlCol && i < pos.brCol && j >= pos.tlRow && j < pos.brRow) {
            var cellIndex = (j-pos.tlRow)*pos.sizeCols + i - pos.tlCol;
            return (cells[cellIndex] || cells[cellIndex].state) ? true : false;
        }

        return false;
    };
};

// Initial patterns (predicate functions for drawing an initial set of live cells)
GOF.initialPatterns = {
    "Default" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,
             1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,
             1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,1,1,1,0,0,
             1,0,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,0,1,0,0,0,0,
             1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,
             1,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,
             0,1,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,1,1,1,1,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,
             1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
             1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,
             1,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1], 25, 25, "center"),
        "sizeCols" : 25,
        "sizeRows" : 25,
        "gridSizeHint": 40
    },
    
    "Gosper Glider Gun" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
             0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
             1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            36, 9, "left-top", 4, 4),
        "sizeCols" : 36,
        "sizeRows" : 9,
        "gridSizeHint" : 50
    },
    
    "Pulsars (oscillators)" : {
        "pfunc" : (function() {
            var cells = [0,0,1,1,1,0,0,0,1,1,1,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         0,0,1,1,1,0,0,0,1,1,1,0,0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,1,1,1,0,0,0,1,1,1,0,0,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         1,0,0,0,0,1,0,1,0,0,0,0,1,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,1,1,1,0,0,0,1,1,1,0,0];
            var sizeCols = 13;
            var sizeRows = 13;
            var centerPulsar = GOF.cellArrayPredicate(cells, sizeCols, sizeRows, "center");
            var pulsars = [GOF.cellArrayPredicate(cells, sizeCols, sizeRows, "center", -13, -13),
                           GOF.cellArrayPredicate(cells, sizeCols, sizeRows, "center", 13, -13),
                           centerPulsar,
                           GOF.cellArrayPredicate(cells, sizeCols, sizeRows, "center", -13, 13),
                           GOF.cellArrayPredicate(cells, sizeCols, sizeRows, "center", 13, 13)];
            return function(i, j, cols, rows, vCols, vRows) {
                if (vCols >= sizeCols*3+1 || vRows >= sizeRows*3+1) {
                    for (var p=0; p<pulsars.length; p++) {
                        if (pulsars[p](i, j, cols, rows, vCols, vRows)) {
                            return true;
                        }
                    }
                    return false;
                } else {
                    return centerPulsar(i,j,cols,rows,vCols,vRows);
                }
            };
        })(),
        "sizeCols": 13,
        "sizeRows": 13,
        "gridSizeHint": 13*3
    },

    "Line of 10 (oscillator)" : {
        "pfunc" : GOF.cellArrayPredicate([1,1,1,1,1,1,1,1,1,1], 10, 1, "center", 0, 0),
        "sizeCols" : 10,
        "sizeRows" : 1,
        "gridSizeHint" : 20
    },
        
    "Random" : {
        "pfunc" : function(i, j, cols, rows, vCols, vRows) {
            return Math.random() < 0.5;
        }
    },

    "2 period oscillators" : {
        "pfunc" : GOF.cellArrayPredicate(
            [1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
             1,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
             0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,
             0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,
             0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,
             0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,
             0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1], 32, 32, "center"),
        "sizeCols":32,
        "sizeRows":32,
        "gridSizeHint": 33
    },

    "R-pentomino" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,1,1,
             1,1,0,
             0,1,0], 3, 3, "center", -10),
        "sizeCols":3,
        "sizeRows":3,
        "gridSizeHint" : 80
    },
    
    "Gosper Centinal (oscillator)" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
             0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
             0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,
             0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,
             0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,
             0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
             1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1], 52, 22),
        "sizeCols": 52,
        "sizeRows": 22,
        "gridSizeHint":54
    },

    "Acorn" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,1,0,0,0,0,0,0,
             0,0,0,1,0,0,0,0,
             1,1,0,0,1,1,1,0], 8, 3, "center"),
        "sizeCols":8,
        "sizeRows":3,
        "gridSizeHint": 50
    },

    "Die hard" : {
        "pfunc": GOF.cellArrayPredicate(
            [0,0,0,0,0,0,1,0,
             1,1,0,0,0,0,0,0,
             0,1,0,0,0,1,1,1], 8, 3, "center"),
        "sizeCols":8,
        "sizeRows":3,
        "gridSizeHint": 25
    },

    "Growth" : {
        "pfunc" : GOF.cellArrayPredicate(
            [1,1,1,0,1,
             1,0,0,0,0,
             0,0,0,1,1,
             0,1,1,0,1,
             1,0,1,0,1], 5, 5, "right-bottom", -20, -20),
        "sizeCols":5,
        "sizeRows":5,
        "gridSizeHint": 110
    },

    "Kok's galaxy (oscillator)" : {
        "pfunc" : GOF.cellArrayPredicate(
            [1,1,0,1,1,1,1,1,1,
             1,1,0,1,1,1,1,1,1,
             1,1,0,0,0,0,0,0,0,
             1,1,0,0,0,0,0,1,1,
             1,1,0,0,0,0,0,1,1,
             1,1,0,0,0,0,0,1,1,
             0,0,0,0,0,0,0,1,1,
             1,1,1,1,1,1,0,1,1,
             1,1,1,1,1,1,0,1,1], 9, 9, "center"),
        "sizeCols":9,
        "sizeRows":9,
        "gridSizeHint" : 14
    },

    "Intersecting circles" : {
        "pfunc" : (function() {
            var predicates = null;
            var grid = null;
            return function(i, j, cols, rows, vCols, vRows) {
                if (!grid ||
                    !(grid.cols == cols && grid.rows == rows && grid.vCols == vCols && grid.vRows == vRows)) {
                    grid = { "cols" : cols, "rows" : rows, "vCols" : vCols, "vRows" : vRows };
                    predicates = new Array();
                    var radius = Math.round(Math.min(vCols, vRows)/4) - 2;
                    predicates.push(GOF.circleCellsPredicate(radius, "center"));
                    predicates.push(GOF.circleCellsPredicate(radius, "center", -radius, -radius));
                    predicates.push(GOF.circleCellsPredicate(radius, "center", -radius, radius));
                    predicates.push(GOF.circleCellsPredicate(radius, "center", radius, -radius));
                    predicates.push(GOF.circleCellsPredicate(radius, "center", radius,  radius));
                }

                for (var p=0; p<predicates.length; p++) {
                    if (predicates[p](i, j, cols, rows, vCols, vRows)) {
                        return true;
                    }
                }
                return false;
            };
        })(),
    },

    "Circles within circles" : {
        "pfunc" : (function() {
            var predicates = null;
            var grid = null;
            return function(i, j, cols, rows, vCols, vRows) {
                if (!grid ||
                    !(grid.cols == cols && grid.rows == rows && grid.vCols == vCols && grid.vRows == vRows)) {
                    grid = { "cols" : cols, "rows" : rows, "vCols" : vCols, "vRows" : vRows };
                    predicates = new Array();

                    var radius = Math.round(Math.min(vCols, vRows)/2) - 2;
                    for (; radius >= 3; radius -= 3) {
                        predicates.push(GOF.circleCellsPredicate(radius, "center"));
                    }
                }

                for (var p=0; p<predicates.length; p++) {
                    if (predicates[p](i, j, cols, rows, vCols, vRows)) {
                        return true;
                    }
                }
                return false;
            };
        })(),
        "gridSizeHint": 68
    },

    "X" : {
        "pfunc" : function(i, j, cols, rows, vCols, vRows) {
            var size = Math.min(vCols,vRows);
            var p = GOF.position(cols, rows, vCols, vRows, size, size, "center", 0, 0);
            if (i >= p.tlCol && i < p.brCol && j >= p.tlRow && j < p.brRow) {
                if (i - p.tlCol == j - p.tlRow) return true;
                if (i - p.tlCol == size - j + p.tlRow - 1) return true;
            }
            return false;
        },
    },

    "Glider attack" : {
        "pfunc" : GOF.cellArrayPredicate(
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,
             0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,
             0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,
             0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,1,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], 51, 54),
        "sizeCols": 51,
        "sizeRows": 54,
        "gridSizeHint" : 55
    },

    "Vertical stripes with noise" : {
        "pfunc" : function(i, j, cols, rows, vCols, vRows) {
            return i % 5 == 1 || Math.random() < 0.05;
        }
    }
};

// Visual themes
GOF.themes = {
    "Default" : {
        "cellShape": "square",
        "cellColorAlive": "#000000",
        "cellColorDead": "#FFFFFF",
        "bgColor": "#CCCCCC",
        "bgColorEnd": "#FFFFFF" 
    },
    
    "Life" : {
        "cellShape": "square",
        "cellColorAlive": "yellow",
        "cellColorDead": "#909090",
        "bgColor": "#BEBEBE",
    },

    "Orange" : {
        "cellShape": "square",
        "cellColorAlive":"orange",
        "cellColorDead":"white",
        "bgColor":"#EEEEEE"
    },

    "Green Dots" :  {
        "cellShape":"circle",
        "cellColorAlive": "#00FF00",
        "cellColorDead":"black",
        "bgColor":"black",
    },
    
    "Red Dots" :  {
        "cellShape":"circle",
        "cellColorAlive": "#FF0000",
        "cellColorDead":"#330000",
        "bgColor":"#000000",
        "bgColorEnd":"#555555"
    },

    "Norway" : {
        "cellShape": "square",
        "cellColorAlive":"#0000FF",
        "cellColorDead":"#FFFFFF",
        "bgColor":"#FF0000"
    },

    "Swedish Meatballs" : {
        "cellShape": "circle",
        "cellColorAlive":"blue",
        "cellColorDead":"yellow",
        "bgColor": "yellow"
    },

    "Default Dots" : {
        "cellShape": "circle",
        "cellColorAlive": "#000000",
        "cellColorDead": "#FFFFFF",
        "bgColor": "#DDDDDD",
        "bgColorEnd": "#FFFFFF" 
    },

};

