// Core Math Engine: GF(4) Arithmetic
// Set elements: 0, 1, 2, 3

// Addition in GF(4) is bitwise XOR
const add = (x, y) => x ^ y;

// Multiplication in GF(4) mapped via polynomial irreducible mod 2: x^2 + x + 1
const mulTable = [
    [0, 0, 0, 0],
    [0, 1, 2, 3],
    [0, 2, 3, 1],
    [0, 3, 1, 2]
];
const mul = (x, y) => mulTable[x][y];

// State variables
let draggingLayerIndex = null;
let layerOffsets = [{x: 0, y: 0}, {x: 0, y: 0}, {x: 0, y: 0}];
let latinSquares = {
    1: [],
    2: [],
    3: []
};

// DOM Elements
const affineGrid = document.getElementById('affine-grid');
const selectSlope = document.getElementById('slope');
const selectIntercept = document.getElementById('intercept');
const eqA = document.getElementById('eq-a');
const eqB = document.getElementById('eq-b');

const ls1Grid = document.querySelector('#ls-1 .matrix-grid');
const ls2Grid = document.querySelector('#ls-2 .matrix-grid');
const ls3Grid = document.querySelector('#ls-3 .matrix-grid');

const selectSq1 = document.getElementById('sq1');
const selectSq2 = document.getElementById('sq2');
const verifyBtn = document.getElementById('verify-btn');
const superMatrix = document.getElementById('super-matrix');
const verifyStatus = document.getElementById('verify-status');
const superpositionDisplay = document.getElementById('superposition-result');

const render3DBtn = document.getElementById('render-3d-btn');
const cubeWrapper = document.getElementById('cube-wrapper');
const cubeScene = document.getElementById('cube-scene');
const tuplePopup = document.getElementById('tuple-popup');
const tupleData = document.getElementById('tuple-data');

const cubeValFilter = document.getElementById('cube-val-filter');

// 3D Scene Rotation State
let isDraggingCube = false;
let previousMousePosition = { x: 0, y: 0 };
let cubeRotation = { x: -30, y: -45 }; // L1 (+Z) at Top-Right, L3 (-Z) at Bottom-Left
let cubeZoom = 1;

// Initialize Application
function init() {
    setupGrid();
    generateAllLatinSquares();
    renderAllLatinSquares();

    // Event Listeners
    selectSlope.addEventListener('change', updateAffinePlane);
    selectIntercept.addEventListener('change', updateAffinePlane);
    verifyBtn.addEventListener('click', verifyOrthogonality);
    render3DBtn.addEventListener('click', render3DCube);
    cubeValFilter.addEventListener('change', render3DCube);
    window.addEventListener('resize', updateAffinePlane);

    // Cube Drag Interactions
    cubeWrapper.addEventListener('mousedown', startCubeDrag);
    cubeWrapper.addEventListener('contextmenu', (e) => e.preventDefault()); // Prevent right-click menu
    cubeWrapper.addEventListener('wheel', handleCubeZoom, { passive: false });
    document.addEventListener('mousemove', dragCube);
    document.addEventListener('mouseup', endCubeDrag);

    // Initial Render
    setTimeout(updateAffinePlane, 50); // wait for DOM rendering
}

// 2D Simulator: Setup 4x4 Grid
function setupGrid() {
    affineGrid.innerHTML = '';

    // Create SVG container for curves //
    const svgNS = "http://www.w3.org/2000/svg";
    const svg = document.createElementNS(svgNS, "svg");
    svg.setAttribute("id", "affine-lines-svg");
    affineGrid.appendChild(svg);

    const labels = ["00", "01", "10", "11"];

    // To make Cartesian (x,y) natural, we render y from 3 down to 0, and x from 0 to 3
    for (let y = 3; y >= 0; y--) {
        for (let x = 0; x <= 3; x++) {
            const point = document.createElement('div');
            point.className = 'grid-point';
            point.dataset.x = x;
            point.dataset.y = y;
            point.id = `pt-${x}-${y}`;

            const dot = document.createElement('div');
            dot.className = 'point-dot';

            const coord = document.createElement('div');
            coord.className = 'point-coord';
            coord.innerText = `(${labels[x]},${labels[y]})`;

            point.appendChild(dot);
            point.appendChild(coord);
            affineGrid.appendChild(point);
        }
    }
}

// Update the highlighted points based on y = ax + b
function updateAffinePlane() {
    const eqDisp = document.querySelector('.equation-display span');
    const aVal = selectSlope.value;
    const bVal = selectIntercept.value;

    const labels = ["00", "01", "10", "11"];
    const aText = aVal === 'all' ? 'all' : (aVal === 'inf' ? '&infin;' : labels[parseInt(aVal, 10)]);
    const bText = bVal === 'all' ? 'all' : labels[parseInt(bVal, 10)];

    if (aVal === 'inf') {
        eqDisp.innerHTML = `x = <span id="eq-b">${bText}</span>`;
    } else {
        eqDisp.innerHTML = `y = <span id="eq-a">${aText}</span>x + <span id="eq-b">${bText}</span>`;
    }

    // Clear all points
    document.querySelectorAll('.point-dot').forEach(dot => dot.classList.remove('active'));
    
    const svg = document.getElementById('affine-lines-svg');
    if (svg) svg.innerHTML = ''; // clear previous paths

    const slopes = aVal === 'all' ? [0, 1, 2, 3, 'inf'] : [aVal === 'inf' ? 'inf' : parseInt(aVal, 10)];
    const intercepts = bVal === 'all' ? [0, 1, 2, 3] : [parseInt(bVal, 10)];

    slopes.forEach(a => {
        intercepts.forEach(b => {
            const activeLinePoints = [];
            
            // Assign distinct colors based on slope for clearer visualization
            const aIndex = a === 'inf' ? 4 : parseInt(a, 10);
            const hueVariations = [210, 45, 120, 280, 0]; // Blue, Yellow, Green, Purple, Red
            const hue = hueVariations[aIndex];
            
            if (a === 'inf') {
                // When displaying "all slopes" for a fixed intercept b, we are constructing 
                // the pencil of lines intersecting at (0, b). Thus, the matching vertical line is x = 0.
                const targetX = (aVal === 'all' && bVal !== 'all') ? 0 : parseInt(b, 10);
                for (let y = 0; y < 4; y++) {
                    const pointDiv = document.getElementById(`pt-${targetX}-${y}`);
                    if (pointDiv) {
                        const dot = pointDiv.querySelector('.point-dot');
                        if (dot) dot.classList.add('active');
                        activeLinePoints.push({ x: targetX, y, elem: pointDiv });
                    }
                }
            } else {
                for (let x = 0; x < 4; x++) {
                    const y = add(mul(a, x), b);
                    const pointDiv = document.getElementById(`pt-${x}-${y}`);
                    if (pointDiv) {
                        const dot = pointDiv.querySelector('.point-dot');
                        if (dot) dot.classList.add('active');
                        activeLinePoints.push({ x, y, elem: pointDiv });
                    }
                }
            }
            
            drawAffineCurve(activeLinePoints, true, hue);
        });
    });
}

function drawAffineCurve(points, append = false, hue = null) {
    const svg = document.getElementById('affine-lines-svg');
    if (!svg) return;

    if (!append) svg.innerHTML = ''; // clear previous paths
    if (points.length < 2) return;

    const gridRect = affineGrid.getBoundingClientRect();

    const coords = points.map(p => {
        const rect = p.elem.getBoundingClientRect();
        return {
            px: rect.left - gridRect.left + rect.width / 2,
            py: rect.top - gridRect.top + rect.height / 2
        };
    });

    // Create an SVG curve path (smooth cubic bezier)
    let d = `M ${coords[0].px},${coords[0].py}`;
    for (let i = 0; i < coords.length - 1; i++) {
        const curr = coords[i];
        const next = coords[i + 1];

        // Control points for a smooth curve
        // We calculate horizontal distance to make the curve look natural
        const dist = Math.abs(next.px - curr.px);
        const cp1x = curr.px + dist * 0.5;
        const cp1y = curr.py;
        const cp2x = next.px - dist * 0.5;
        const cp2y = next.py;

        d += ` C ${cp1x},${cp1y} ${cp2x},${cp2y} ${next.px},${next.py}`;
    }

    const pathNS = "http://www.w3.org/2000/svg";
    const path = document.createElementNS(pathNS, "path");
    path.setAttribute("d", d);
    path.setAttribute("class", "affine-curve");

    if (hue !== null) {
        path.style.stroke = `hsl(${hue}, 80%, 65%)`;
        path.style.filter = `drop-shadow(0 0 8px hsl(${hue}, 80%, 50%))`;
    }

    svg.appendChild(path);
}

// Generate Latin Squares
// For a fixed slope a, the value at (x, y) is the intercept b
// b = y - ax => b = add(y, mul(a, x))
function generateLatinSquare(slope) {
    const matrix = [];
    // Typically matrices are drawn row by row. 
    // To map Cartesian to Matrix, a typical way is:
    // Rows top-to-bottom could be mapped to y=0 to 3, columns to x=0 to 3.
    // We will use standard row=y, col=x where y is 0 to 3 from top to bottom.
    // wait, our math coordinates had y=0..3. Let's just iterate y from 3 down to 0 to match visual grid, 
    // or just iterate row from 0 to 3 sequentially. Mathematical definition doesn't enforce visual row order.
    // We'll iterate row=0..3, col=0..3 => Let's map x = col, y = row
    for (let r = 0; r < 4; r++) {
        let rowData = [];
        for (let c = 0; c < 4; c++) {
            const x = c;
            const y = r;
            const intercept_b = add(y, mul(slope, x));
            rowData.push(intercept_b);
        }
        matrix.push(rowData);
    }
    return matrix;
}

function generateAllLatinSquares() {
    latinSquares[1] = generateLatinSquare(1);
    latinSquares[2] = generateLatinSquare(2);
    latinSquares[3] = generateLatinSquare(3);
}

function renderMatrix(gridElement, matrix) {
    gridElement.innerHTML = '';
    for (let r = 0; r < 4; r++) {
        for (let c = 0; c < 4; c++) {
            const cell = document.createElement('div');
            cell.className = 'matrix-cell';

            // Assign a color class or direct style based on value for visual appeal
            const val = matrix[r][c];
            cell.innerText = val;

            // Optional: Give different hues to different numbers
            const hue = [0, 210, 45, 120][val];
            cell.style.color = `hsl(${hue}, 80%, 75%)`;

            gridElement.appendChild(cell);
        }
    }
}

function renderAllLatinSquares() {
    renderMatrix(ls1Grid, latinSquares[1]);
    renderMatrix(ls2Grid, latinSquares[2]);
    renderMatrix(ls3Grid, latinSquares[3]);
}

// Superimposition Verification
function verifyOrthogonality() {
    const s1Id = parseInt(selectSq1.value, 10);
    const s2Id = parseInt(selectSq2.value, 10);

    if (s1Id === s2Id) {
        showVerificationResult(false, `Cannot superimpose a square with itself! Select different squares.`);
        return;
    }

    const mat1 = latinSquares[s1Id];
    const mat2 = latinSquares[s2Id];

    superMatrix.innerHTML = '';
    const tupleSet = new Set();

    for (let r = 0; r < 4; r++) {
        for (let c = 0; c < 4; c++) {
            const val1 = mat1[r][c];
            const val2 = mat2[r][c];
            const tuplePair = `${val1}${val2}`;

            tupleSet.add(tuplePair);

            const cell = document.createElement('div');
            cell.className = 'matrix-cell';
            cell.innerText = tuplePair;

            // Subtle color mix 
            const hue1 = [0, 210, 45, 120][val1];
            const hue2 = [0, 210, 45, 120][val2];
            cell.style.background = `linear-gradient(135deg, hsla(${hue1}, 80%, 30%, 0.3), hsla(${hue2}, 80%, 30%, 0.3))`;

            superMatrix.appendChild(cell);
        }
    }

    superpositionDisplay.classList.remove('hidden');

    // Verifying if Set size is exactly 16
    if (tupleSet.size === 16) {
        showVerificationResult(true, `驗證成功：互為正交 (Mutually Orthogonal)<br><span style="font-size:0.9rem;font-weight:normal;opacity:0.8;margin-top:5px;display:block">Found exactly ${tupleSet.size} unique pairs.</span>`);
    } else {
        showVerificationResult(false, `驗證失敗：非正交 (Not Orthogonal)<br><span style="font-size:0.9rem;font-weight:normal;opacity:0.8;margin-top:5px;display:block">Found only ${tupleSet.size} unique pairs.</span>`);
    }
}

function showVerificationResult(success, htmlMessage) {
    verifyStatus.innerHTML = '';
    verifyStatus.className = 'verification-status'; // Reset

    const icon = document.createElement('span');
    icon.className = 'status-icon';
    icon.innerHTML = success ? '&check;' : '&times;';

    const msg = document.createElement('div');
    msg.innerHTML = htmlMessage;

    if (success) {
        verifyStatus.classList.add('status-success');
    } else {
        verifyStatus.classList.add('status-error');
    }

    verifyStatus.appendChild(icon);
    verifyStatus.appendChild(msg);
}

// 3D Overlaid Planes Representation
function render3DCube() {
    cubeWrapper.classList.remove('hidden');
    cubeScene.innerHTML = '';
    tuplePopup.classList.add('hidden'); // Hide popup on re-render

    // Make sure transform respects initial state
    cubeScene.style.transform = `scale(${cubeZoom}) rotateX(${cubeRotation.x}deg) rotateY(${cubeRotation.y}deg)`;

    // We stack L1, L2, L3 as separated planes in the Z-axis.
    const layers = [1, 2, 3];
    const spacing = 120; // Z-spacing between layers
    const offsetZ = 1.0 * spacing;
    const labels = ["00", "01", "10", "11"];

    // Read filter state
    const filterValStr = cubeValFilter.value;
    const targetVal = filterValStr === 'all' ? 'all' : parseInt(filterValStr, 10);

    // Reset layer offsets
    layerOffsets = [{x: 0, y: 0}, {x: 0, y: 0}, {x: 0, y: 0}];

    layers.forEach((slope, zIndex) => {
        const matrix = latinSquares[slope];

        // Create a layer plane
        const layerPlane = document.createElement('div');
        layerPlane.className = 'layer-plane';

        // Position the plane in 3D space along the Z-axis
        const posZ = offsetZ - (zIndex * spacing);
        layerPlane.style.transform = `translate3d(0px, 0px, ${posZ}px)`;
        layerPlane.dataset.initZ = posZ;

        // Background label for the plane
        const layerLabel = document.createElement('div');
        layerLabel.className = 'layer-label';
        layerLabel.innerText = `L${slope}`;
        layerPlane.appendChild(layerLabel);

        for (let r = 0; r < 4; r++) {
            for (let c = 0; c < 4; c++) {
                const val = matrix[r][c];

                const v1 = latinSquares[1][r][c];
                const v2 = latinSquares[2][r][c];
                const v3 = latinSquares[3][r][c];
                const tupleString = `(${labels[v1]}, ${labels[v2]}, ${labels[v3]})`;

                // Create a 2D cell in the plane
                const cell = document.createElement('div');
                cell.className = 'layer-cell';

                // Color coding based on value: dark center, light edges
                const hue = [0, 210, 45, 120][val];
                cell.style.background = `radial-gradient(circle at center, hsla(${hue}, 80%, 20%, 0.95) 0%, hsla(${hue}, 80%, 55%, 0.65) 100%)`;
                cell.innerText = labels[val];
                
                // Hide blocks that do not belong to the selected value
                if (targetVal === 'all' || val === targetVal) {
                    cell.style.opacity = '1';
                } else {
                    cell.style.opacity = '0.05'; // 'Ghost' mode for unselected geometric space
                }

                // Click interaction
                cell.addEventListener('click', (e) => {
                    e.stopPropagation();
                    tuplePopup.classList.remove('hidden');
                    tupleData.innerText = tupleString;

                    // Add a tiny animation feedback to the clicked cell
                    cell.style.transform = `translateZ(15px) scale(1.05)`;
                    setTimeout(() => {
                        cell.style.transform = ``;
                    }, 300);
                });

                layerPlane.appendChild(cell);
            }
        }

        cubeScene.appendChild(layerPlane);
    });
}

// Cube Rotation Logic
function startCubeDrag(e) {
    if (e.button === 0) {
        // Left click: Layer dragging
        const plane = e.target.closest('.layer-plane');
        if (plane) {
            draggingLayerIndex = Array.from(cubeScene.children).indexOf(plane);
            isDraggingCube = false;
        } else {
            // Clicking empty space with left mouse button does nothing
            draggingLayerIndex = null;
            isDraggingCube = false;
        }
    } else if (e.button === 2) {
        // Right click: Scene rotation
        isDraggingCube = true;
        draggingLayerIndex = null;
    } else {
        return; // Ignore middle click or other buttons
    }
    
    previousMousePosition = { x: e.offsetX, y: e.offsetY };
}

function dragCube(e) {
    const deltaMove = {
        x: e.movementX || e.mozMovementX || e.webkitMovementX || 0,
        y: e.movementY || e.mozMovementY || e.webkitMovementY || 0
    };

    if (draggingLayerIndex !== null) {
        layerOffsets[draggingLayerIndex].x += deltaMove.x;
        layerOffsets[draggingLayerIndex].y += deltaMove.y;
        
        const targetPlane = cubeScene.children[draggingLayerIndex];
        const initZ = targetPlane.dataset.initZ;
        targetPlane.style.transform = `translate3d(${layerOffsets[draggingLayerIndex].x}px, ${layerOffsets[draggingLayerIndex].y}px, ${initZ}px)`;
    } else if (isDraggingCube) {
        // Calculate rotation: dragging X rotates around Y axis, dragging Y rotates around X axis
        cubeRotation.y += deltaMove.x * 0.5;
        cubeRotation.x -= deltaMove.y * 0.5;

        // Clamp X rotation to avoid flipping upside down entirely for better orientation feeling
        cubeRotation.x = Math.max(-80, Math.min(80, cubeRotation.x));

        cubeScene.style.transform = `scale(${cubeZoom}) rotateX(${cubeRotation.x}deg) rotateY(${cubeRotation.y}deg)`;
    }
}

function handleCubeZoom(e) {
    if (!cubeWrapper.classList.contains('hidden')) {
        e.preventDefault(); // Prevent page scroll when zooming inside 3D viewer
        const zoomSpeed = 0.05;
        const zoomDelta = e.deltaY > 0 ? -zoomSpeed : zoomSpeed;
        cubeZoom += zoomDelta;
        
        // Clamp scale to reasonable bounds
        cubeZoom = Math.max(0.3, Math.min(3.0, cubeZoom));
        
        cubeScene.style.transform = `scale(${cubeZoom}) rotateX(${cubeRotation.x}deg) rotateY(${cubeRotation.y}deg)`;
    }
}

function endCubeDrag(e) {
    isDraggingCube = false;
    draggingLayerIndex = null;
}

// Start app
document.addEventListener('DOMContentLoaded', init);
