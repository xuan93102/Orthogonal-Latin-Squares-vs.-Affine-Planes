let THREE, OrbitControls;

function formatAlphaLabel(coeffs) {
    const superscripts = ['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'];
    const hasNonZero = coeffs.some(c => c !== 0);
    if (!hasNonZero) return '0';

    const parts = [];
    for (let degree = coeffs.length - 1; degree >= 0; degree--) {
        const coeff = coeffs[degree];
        if (coeff === 0) continue;
        if (degree === 0) {
            parts.push(`${coeff}`);
            continue;
        }
        const power = degree === 1 ? 'α' : `α${String(degree).split('').map(d => superscripts[parseInt(d, 10)]).join('')}`;
        if (coeff === 1) {
            parts.push(power);
        } else {
            parts.push(`${coeff}${power}`);
        }
    }
    return parts.join(' + ');
}

function createGFPrime(p) {
    return {
        q: p,
        p,
        add: (a, b) => (a + b) % p,
        mul: (a, b) => (a * b) % p,
        label: (i) => i.toString().padStart(2, '0'),
        type: 'prime',
        desc: `質數體 GF(${p})\n加法: (a+b) mod ${p}\n乘法: (a×b) mod ${p}`
    };
}

function createGF2k(q, k) {
    const polys = { 1: 3, 2: 7, 3: 11, 4: 19 };
    const poly = polys[k];

    function mul(a, b) {
        let res = 0;
        let tempA = a;
        for (let i = 0; i < k; i++) {
            if ((b >> i) & 1) res ^= tempA;
            let hi = tempA & (1 << (k - 1));
            tempA <<= 1;
            if (hi) tempA ^= poly;
        }
        return res;
    }

    const intToBinaryCoeffs = (value) => {
        const coeffs = new Array(k).fill(0);
        for (let bit = 0; bit < k; bit++) {
            coeffs[bit] = (value >> bit) & 1;
        }
        return coeffs;
    };

    return {
        q,
        p: 2,
        add: (a, b) => a ^ b,
        mul: mul,
        label: (i) => i.toString(2).padStart(k, '0'),
        type: 'binary',
        desc: `擴張體 GF(2^${k}) = GF(${q})\n加法: 多項式 XOR (逐位元 XOR)\n乘法: 模不可約多項式算術`
    };
}

function createGFpk(p, k) {
    if (k === 1) return createGFPrime(p);
    if (p === 2) return createGF2k(Math.pow(2, k), k);

    const irreduciblePolys = {
        '3,2': [1, 0, 1], // x² + 1 over GF(3)
    };
    const poly = irreduciblePolys[`${p},${k}`];
    if (!poly) return null;

    const q = Math.pow(p, k);
    const mod = n => ((n % p) + p) % p;

    const intToCoeffs = (value) => {
        const coeffs = new Array(k).fill(0);
        let x = value;
        for (let i = 0; i < k; i++) {
            coeffs[i] = x % p;
            x = Math.floor(x / p);
        }
        return coeffs;
    };

    const coeffsToInt = (coeffs) => coeffs.reduce((acc, c, idx) => acc + c * Math.pow(p, idx), 0);

    const add = (a, b) => {
        const ca = intToCoeffs(a);
        const cb = intToCoeffs(b);
        return coeffsToInt(ca.map((value, index) => mod(value + cb[index])));
    };

    const mul = (a, b) => {
        const ca = intToCoeffs(a);
        const cb = intToCoeffs(b);
        const product = new Array(2 * k - 1).fill(0);

        for (let i = 0; i < k; i++) {
            for (let j = 0; j < k; j++) {
                product[i + j] = mod(product[i + j] + ca[i] * cb[j]);
            }
        }

        for (let deg = product.length - 1; deg >= k; deg--) {
            const coefficient = product[deg];
            if (coefficient === 0) continue;
            for (let i = 0; i <= k; i++) {
                product[deg - (k - i)] = mod(product[deg - (k - i)] - coefficient * poly[i]);
            }
            product[deg] = 0;
        }

        return coeffsToInt(product.slice(0, k).map(mod));
    };

    const label = (value) => formatAlphaLabel(intToCoeffs(value));

    const polyDescription = poly
        .map((coef, idx) => {
            const degree = poly.length - 1 - idx;
            if (coef === 0) return null;
            if (degree === 0) return `${coef}`;
            if (coef === 1) return degree === 1 ? 'x' : `x^${degree}`;
            return degree === 1 ? `${coef}x` : `${coef}x^${degree}`;
        })
        .filter(Boolean)
        .join(' + ');

    return {
        q,
        p,
        add,
        mul,
        label,
        type: 'extension',
        desc: `擴張體 GF(${p}^${k}) = GF(${q})\n加法: 多項式係數逐項 mod ${p}\n乘法: 模不可約多項式 ${polyDescription}`
    };
}

function createField(q, p = null, k = null) {
    if (p !== null && k !== null) {
        const field = createGFpk(p, k);
        if (field) return field;
    }

    if (q === 2) return createGF2k(2, 1);
    if (q === 4) return createGF2k(4, 2);
    if (q === 8) return createGF2k(8, 3);
    if (q === 16) return createGF2k(16, 4);

    const primes = [2, 3, 5, 7, 11, 13];
    if (primes.includes(q)) return createGFPrime(q);

    return null;
}

const VALID_Q = [2, 3, 4, 5, 7, 8, 11, 13, 16];

const PALETTE = [
    '#5fa2ce', '#b48ead', '#d08770', '#a3be8c',
    '#ebcb8b', '#bf616a', '#88c0d0', '#81a1c1',
    '#f4a261', '#e76f51', '#2a9d8f', '#e9c46a',
    '#264653', '#6a4c93', '#c77dff', '#48cae4',
    '#f77f00', '#d62828', '#023e8a', '#80b918',
    '#bc6c25', '#606c38', '#ddb892', '#9b5de5',
    '#00bbf9',
];

let currentField = null;
let currentA = 1;
let currentB = 1;
let show3D = false; // true if GF(2^k)
let showNumber = false;

let svgEl, latinSquareEl, eqDisplay2D, eqDisplay3D, aSlider, aSliderVal, fieldDesc;
let bSlider, bSliderVal;
let kInput, pkPInput, runKBtn;
let controls3d, showNumberCheck, view2d, view3d;
let axisXSelect, axisYSelect, axisZSelect;
let selectedAxisX = 'all', selectedAxisY = 'all', selectedAxisZ = 'all';

let scene3dLines, scene3dLatin; 

class Visualizer3D {
    constructor(container) {
        this.container = container;
        this.init();
    }

    init() {
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(45, this.container.clientWidth / 500, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        this.renderer.setSize(this.container.clientWidth, 500);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.camera.position.set(12, 12, 12);
        this.controls.update();

        const light = new THREE.DirectionalLight(0xffffff, 3);
        light.position.set(5, 10, 7.5);
        this.scene.add(light);
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.5));

        this.objects = new THREE.Group();
        this.scene.add(this.objects);

        this.animate = this.animate.bind(this);
        this.animate();

        window.addEventListener('resize', () => {
            if (this.container.clientWidth === 0) return;
            this.camera.aspect = this.container.clientWidth / 500;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(this.container.clientWidth, 500);
        });
    }

    animate() {
        requestAnimationFrame(this.animate);
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    clear() {
        while (this.objects.children.length > 0) {
            const obj = this.objects.children[0];
            if (obj.geometry) obj.geometry.dispose();
            if (obj.material) {
                if (Array.isArray(obj.material)) obj.material.forEach(m => m.dispose());
                else obj.material.dispose();
            }
            this.objects.remove(obj);
        }
    }

    render(mode, field, showText, aVal, bVal, filters = { x: 'all', y: 'all', z: 'all' }) {
        this.clear();
        const q = field.q;
        const spacing = 1.6;
        const offset = (q - 1) * spacing / 2;

        this.drawAxes(q, spacing, offset, field, mode);

        const xFilter = filters.x === 'all' ? null : Number(filters.x);
        const yFilter = filters.y === 'all' ? null : Number(filters.y);
        const zFilter = filters.z === 'all' ? null : Number(filters.z);

        if (mode === '3d_latin') {
            const group = new THREE.Group();
            for (let c = 0; c < q; c++) {
                const color = PALETTE[c % PALETTE.length];
                const mat = new THREE.MeshPhongMaterial({
                    color: new THREE.Color(color), transparent: true, opacity: 0.65
                });

                for (let x = 0; x < q; x++) {
                    if (xFilter !== null && x !== xFilter) continue;
                    for (let y = 0; y < q; y++) {
                        if (yFilter !== null && y !== yFilter) continue;
                        const ax = field.mul(aVal, x);
                        const by = field.mul(bVal, y);
                        const z = field.add(field.add(ax, by), c);
                        if (zFilter !== null && z !== zFilter) continue;

                        const boxGeom = new THREE.BoxGeometry(0.9, 0.9, 0.08);
                        const mesh = new THREE.Mesh(boxGeom, mat);
                        mesh.position.set(x * spacing - offset, y * spacing - offset, z * spacing - offset);
                        group.add(mesh);

                        if (showText) {
                            const sprite = this.makeTextSprite(field.label(c), "white", 44);
                            sprite.position.set(x * spacing - offset, y * spacing - offset + 0.5, z * spacing - offset + 0.1);
                            sprite.scale.set(0.75, 0.75, 1);
                            group.add(sprite);
                        }
                    }
                }
            }
            this.objects.add(group);
        } else if (mode === '3d_lines') {
            for (let c = 0; c < q; c++) {
                const color = PALETTE[c % PALETTE.length];
                const group = new THREE.Group();
                const material = new THREE.LineBasicMaterial({
                    color: color, transparent: true, opacity: 0.9, linewidth: 3
                });
                const dotMat = new THREE.MeshBasicMaterial({ color: color });
                const dotGeom = new THREE.SphereGeometry(0.12, 12, 12);

                const drawnDots = new Set();
                const drawDot = (x, y, z) => {
                    const id = `${x},${y},${z}`;
                    if (!drawnDots.has(id)) {
                        drawnDots.add(id);
                        const dot = new THREE.Mesh(dotGeom, dotMat);
                        dot.position.set(x * spacing - offset, y * spacing - offset, z * spacing - offset);
                        group.add(dot);
                    }
                };

                for (let y = 0; y < q; y++) {
                    if (yFilter !== null && y !== yFilter) continue;
                    const points = [];
                    for (let x = 0; x < q; x++) {
                        if (xFilter !== null && x !== xFilter) continue;
                        const ax = field.mul(aVal, x);
                        const by = field.mul(bVal, y);
                        const z = field.add(field.add(ax, by), c);
                        if (zFilter !== null && z !== zFilter) continue;
                        points.push(new THREE.Vector3(x * spacing - offset, y * spacing - offset, z * spacing - offset));
                        drawDot(x, y, z);
                    }
                    if (points.length <= 1) continue;
                    if (aVal === 0) {
                        group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
                    } else {
                        group.add(new THREE.Line(this.makeCurveGeometry(points, points.length, 'x'), material));
                    }
                }

                for (let x = 0; x < q; x++) {
                    if (xFilter !== null && x !== xFilter) continue;
                    const points = [];
                    for (let y = 0; y < q; y++) {
                        if (yFilter !== null && y !== yFilter) continue;
                        const ax = field.mul(aVal, x);
                        const by = field.mul(bVal, y);
                        const z = field.add(field.add(ax, by), c);
                        if (zFilter !== null && z !== zFilter) continue;
                        points.push(new THREE.Vector3(x * spacing - offset, y * spacing - offset, z * spacing - offset));
                    }
                    if (points.length <= 1) continue;
                    if (bVal === 0) {
                        group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
                    } else {
                        group.add(new THREE.Line(this.makeCurveGeometry(points, points.length, 'y'), material));
                    }
                }

                this.objects.add(group);
            }
        }
    }

    makeCurveGeometry(points, q, direction) {
        const curvePath = new THREE.CurvePath();
        for (let k = 0; k < points.length - 1; k++) {
            const p0 = points[k];
            const p1 = points[k + 1];
            const cpd = direction === 'x' ? (p1.x - p0.x) * 0.42 : (p1.y - p0.y) * 0.42;
            const cp1 = new THREE.Vector3(
                direction === 'x' ? p0.x + cpd : p0.x, 
                direction === 'y' ? p0.y + cpd : p0.y, 
                p0.z
            );
            const cp2 = new THREE.Vector3(
                direction === 'x' ? p1.x - cpd : p1.x, 
                direction === 'y' ? p1.y - cpd : p1.y, 
                p1.z
            );
            curvePath.add(new THREE.CubicBezierCurve3(p0, cp1, cp2, p1));
        }
        return new THREE.BufferGeometry().setFromPoints(curvePath.getPoints(q * 8));
    }

    drawAxes(q, spacing, offset, field, mode) {
        const len = q * spacing;
        const o = -offset - 1.8;
        const origin = new THREE.Vector3(o, o, o);

        this.objects.add(new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), origin, len + 3, 0xff4444));
        this.objects.add(new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), origin, len + 3, 0x44ff44));
        this.objects.add(new THREE.ArrowHelper(new THREE.Vector3(0, 0, 1), origin, len + 3, 0x4444ff));

        const lxTitle = this.makeTextSprite('X (x)', "#ff8888", 60);
        lxTitle.position.set(o + len + 4, o, o);
        lxTitle.scale.set(2, 2, 1);
        this.objects.add(lxTitle);

        const lyTitle = this.makeTextSprite('Y (y)', "#88ff88", 60);
        lyTitle.position.set(o, o + len + 4, o);
        lyTitle.scale.set(2, 2, 1);
        this.objects.add(lyTitle);

        const lzTitle = this.makeTextSprite('Z (z)', "#8888ff", 60);
        lzTitle.position.set(o, o, o + len + 4);
        lzTitle.scale.set(2, 2, 1);
        this.objects.add(lzTitle);

        for (let i = 0; i < q; i++) {
            const lbl = field.label(i);
            const size = 1.1;

            const lx = this.makeTextSprite(lbl, "#ffaaaa", 40);
            lx.position.set(i * spacing - offset, o - 0.8, o);
            lx.scale.set(size, size, 1);
            this.objects.add(lx);

            const ly = this.makeTextSprite(lbl, "#aaffaa", 40);
            ly.position.set(o - 0.8, i * spacing - offset, o);
            ly.scale.set(size, size, 1);
            this.objects.add(ly);

            const lz = this.makeTextSprite(lbl, "#aaaaff", 40);
            lz.position.set(o, o - 0.8, i * spacing - offset);
            lz.scale.set(size, size, 1);
            this.objects.add(lz);
        }
    }

    makeTextSprite(message, color = "white", fontSize = 42) {
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        canvas.width = 128;
        canvas.height = 128;
        ctx.fillStyle = color;
        ctx.font = `bold ${fontSize}px sans-serif`;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(message, 64, 64);

        const texture = new THREE.CanvasTexture(canvas);
        const spriteMaterial = new THREE.SpriteMaterial({ map: texture, depthTest: false });
        const sprite = new THREE.Sprite(spriteMaterial);
        sprite.renderOrder = 999;
        return sprite;
    }
}

function init() {
    svgEl = document.getElementById('draw-canvas');
    latinSquareEl = document.getElementById('latin-square');
    eqDisplay2D = document.getElementById('eq-display-2d');
    eqDisplay3D = document.getElementById('eq-display-3d');
    aSlider = document.getElementById('a-slider');
    aSliderVal = document.getElementById('a-slider-val');
    bSlider = document.getElementById('b-slider');
    bSliderVal = document.getElementById('b-slider-val');
    fieldDesc = document.getElementById('field-desc');

    kInput = document.getElementById('k-input');
    pkPInput = document.getElementById('pk-p-input');
    runKBtn = document.getElementById('run-k');

    controls3d = document.getElementById('controls-3d');
    showNumberCheck = document.getElementById('show-number-check');
    axisXSelect = document.getElementById('axis-x-select');
    axisYSelect = document.getElementById('axis-y-select');
    axisZSelect = document.getElementById('axis-z-select');
    view2d = document.getElementById('view-2d');
    view3d = document.getElementById('view-3d');

    scene3dLines = new Visualizer3D(document.getElementById('viz-3d-lines'));
    scene3dLatin = new Visualizer3D(document.getElementById('viz-3d-latin'));

    runKBtn.onclick = () => {
        const p = parseInt(pkPInput.value);
        const k = parseInt(kInput.value);
        const maxK = getMaxKForPrime(p);
        const q = Math.pow(p, k);
        if (!isPrime(p)) {
            alert("請選擇有效的質數 p");
            return;
        }
        if (k < 1 || k > maxK) {
            alert(`k 必須在 1 到 ${maxK} 之間，且 q = p^k 不可超過 16`);
            return;
        }
        switchField(q, p, k);
    };

    pkPInput.addEventListener('change', () => updateKLimit());
    updateKLimit();

    showNumberCheck.onchange = e => {
        showNumber = e.target.checked;
        if (show3D) render();
    };

    axisXSelect.onchange = e => {
        selectedAxisX = e.target.value;
        render();
    };
    axisYSelect.onchange = e => {
        selectedAxisY = e.target.value;
        render();
    };
    axisZSelect.onchange = e => {
        selectedAxisZ = e.target.value;
        render();
    };

    aSlider.addEventListener('input', e => {
        const v = parseInt(e.target.value);
        setA(v === -1 ? 'inf' : v);
    });

    bSlider.addEventListener('input', e => {
        const v = parseInt(e.target.value);
        setB(v);
    });

    switchField(4);
}

function isPrime(n) {
    if (n < 2) return false;
    for (let i = 2; i <= Math.sqrt(n); i++) if (n % i === 0) return false;
    return true;
}

function getMaxKForPrime(p) {
    return Math.max(1, Math.floor(Math.log(16) / Math.log(p)));
}

function updateKLimit() {
    if (!pkPInput) return;
    const p = parseInt(pkPInput.value);
    const maxK = getMaxKForPrime(p);
    kInput.max = maxK;
    if (parseInt(kInput.value) > maxK) kInput.value = maxK;
}

function updateDisplayMode() {
    if (show3D) {
        view2d.style.display = 'flex';
        view3d.style.display = 'flex';
        document.getElementById('b-slider-group').style.display = 'flex';
        document.getElementById('a-slider-label-text').innerText = '選擇 X 係數 (a)';
        if (currentA === -1) setA(0);
        aSlider.min = 0;
        eqDisplay3D.style.display = 'block';
        setTimeout(() => window.dispatchEvent(new Event('resize')), 10);
    } else {
        view2d.style.display = 'flex';
        view3d.style.display = 'none';
        document.getElementById('b-slider-group').style.display = 'none';
        document.getElementById('a-slider-label-text').innerText = '選擇斜率 (a 值)';
        aSlider.min = -1;
        eqDisplay3D.style.display = 'none';
    }
    render();
}

function switchField(q) {
    currentField = createField(q);
    if (!currentField) return;

    currentA = 1;
    currentB = 1;
    selectedAxisX = 'all';
    selectedAxisY = 'all';
    selectedAxisZ = 'all';

    controls3d.style.display = 'flex';
    show3D = true;

    aSlider.min = 0;
    aSlider.max = q - 1;
    aSlider.value = currentA;

    updateAxisSelectors();

    bSlider.min = 0;
    bSlider.max = q - 1;
    bSlider.value = currentB;

    fieldDesc.innerText = currentField.desc;
    updateSliderLabelA();
    updateSliderLabelB();
    updateDisplayMode();
}

function setA(val) {
    currentA = (val === 'inf') ? -1 : parseInt(val);
    aSlider.value = currentA;
    updateSliderLabelA();
    render();
}

function setB(val) {
    currentB = parseInt(val);
    bSlider.value = currentB;
    updateSliderLabelB();
    render();
}

function updateAxisSelectors() {
    if (!currentField || !axisXSelect || !axisYSelect || !axisZSelect) return;
    const q = currentField.q;
    const isBinary = currentField.type === 'binary';
    const bitWidth = Math.log2(q);

    const makeOption = (value, label) => {
        const opt = document.createElement('option');
        opt.value = value;
        opt.textContent = label;
        return opt;
    };

    [axisXSelect, axisYSelect, axisZSelect].forEach(select => {
        select.innerHTML = '';
        select.appendChild(makeOption('all', 'All'));
        for (let i = 0; i < q; i++) {
            const label = isBinary
                ? i.toString(2).padStart(bitWidth, '0')
                : currentField.label(i);
            select.appendChild(makeOption(String(i), label));
        }
    });

    axisXSelect.value = selectedAxisX;
    axisYSelect.value = selectedAxisY;
    axisZSelect.value = selectedAxisZ;
}

function updateSliderLabelA() {
    aSliderVal.textContent = (currentA === -1) ? 'a = ∞ (垂直線)' : `a = ${currentField.label(currentA)}`;
}

function updateSliderLabelB() {
    bSliderVal.textContent = `b = ${currentField.label(currentB)}`;
}

function render() {
    if (show3D) {
        eqDisplay3D.innerHTML = `3D: z = ${currentField.label(currentA)}·x + ${currentField.label(currentB)}·y + c`;
        const axisFilters = {
            x: selectedAxisX,
            y: selectedAxisY,
            z: selectedAxisZ
        };
        scene3dLines.render('3d_lines', currentField, showNumber, currentA, currentB, axisFilters);
        scene3dLatin.render('3d_latin', currentField, showNumber, currentA, currentB, axisFilters);
    }
    render2D();
    renderLatinSquare(currentField, currentField.q);
}

function render2D() {
    const field = currentField;
    const q = field.q;

    eqDisplay2D.textContent = (currentA === -1) ? `2D: x = c` : `2D: y = ${field.label(currentA)} · x + c`;

    const maxPlot = 380;
    const numGaps = q - 1;
    const gridSize = numGaps > 0 ? Math.round(maxPlot / numGaps) : maxPlot;

    const fontSize = Math.max(8, Math.min(14, gridSize * 0.28));
    const dotR = Math.max(2, Math.min(7, gridSize * 0.10));
    const strokeW = Math.max(1, Math.min(4, gridSize * 0.06));

    const sideMargin = Math.max(80, fontSize * 6);
    const topMargin = 40;
    const bottomMargin = Math.max(50, fontSize * 4);

    const svgW = sideMargin + (numGaps * gridSize) + sideMargin;
    const svgH = topMargin + (numGaps * gridSize) + bottomMargin;

    const offX = sideMargin;
    const offY = topMargin;

    svgEl.setAttribute('viewBox', `0 0 ${svgW} ${svgH}`);
    svgEl.innerHTML = '';

    const px = xi => offX + xi * gridSize;
    const py = yi => offY + (numGaps - yi) * gridSize;

    for (let i = 0; i < q; i++) {
        const lbl = field.label(i);
        // X 軸
        svgEl.appendChild(makeText(lbl, px(i), offY + numGaps * gridSize + fontSize + 8, 'axis-label', fontSize));
        // Y 軸
        const ty = makeText(lbl, offX - 12, py(i) + fontSize * 0.38, 'axis-label', fontSize);
        ty.setAttribute('text-anchor', 'end');
        svgEl.appendChild(ty);
    }

    for (let b = 0; b < q; b++) {
        let pts = [];
        if (currentA === -1) {
            for (let y = 0; y < q; y++) pts.push({ x: b, y });
        } else {
            for (let x = 0; x < q; x++) {
                const y = field.add(field.mul(currentA, x), b);
                pts.push({ x, y });
            }
        }
        pts.sort((p1, p2) => p1.x - p2.x);

        const color = PALETTE[b % PALETTE.length];
        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        g.setAttribute('class', 'line-group');

        let isStraight = (currentA === -1 || currentA === 0);

        let d = `M ${px(pts[0].x)} ${py(pts[0].y)}`;
        if (isStraight) {
            for (let k = 1; k < pts.length; k++) d += ` L ${px(pts[k].x)} ${py(pts[k].y)}`;
        } else {
            for (let k = 0; k < pts.length - 1; k++) {
                const x0 = px(pts[k].x), y0 = py(pts[k].y);
                const x1 = px(pts[k + 1].x), y1 = py(pts[k + 1].y);
                const cpd = (x1 - x0) * 0.42;
                d += ` C ${x0 + cpd} ${y0}, ${x1 - cpd} ${y1}, ${x1} ${y1}`;
            }
        }

        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path.setAttribute('d', d);
        path.setAttribute('class', 'line-path');
        path.setAttribute('stroke', color);
        path.setAttribute('stroke-width', strokeW);
        g.appendChild(path);

        // 懸停標籤
        const last = pts[pts.length - 1];
        const lblStr = (currentA === -1) ? `x=${field.label(b)}` : `b=${field.label(b)}`;
        const lbl = makeText(lblStr, px(last.x) + 7, py(last.y) + fontSize * 0.3, 'equation-label', Math.max(8, fontSize - 1));
        lbl.setAttribute('fill', color);
        lbl.setAttribute('text-anchor', 'start');
        g.appendChild(lbl);

        svgEl.appendChild(g);
    }

    for (let y = 0; y < q; y++) {
        for (let x = 0; x < q; x++) {
            const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            c.setAttribute('cx', px(x));
            c.setAttribute('cy', py(y));
            c.setAttribute('r', dotR);
            c.setAttribute('class', 'grid-point');
            svgEl.appendChild(c);
        }
    }
}

function renderLatinSquare(field, q) {
    latinSquareEl.innerHTML = '';
    latinSquareEl.style.gridTemplateColumns = `repeat(${q}, 1fr)`;

    const cellSize = Math.max(20, Math.min(60, Math.floor(280 / q)));
    const cellFont = Math.max(7, Math.min(16, cellSize * 0.4));
    latinSquareEl.style.gap = Math.max(1, Math.min(5, Math.floor(cellSize * 0.1))) + 'px';

    if (currentA === -1 || currentA === 0) {
        for (let i = 0; i < q * q; i++) {
            latinSquareEl.appendChild(makeCell('−', '#444', cellSize, cellFont));
        }
    } else {
        for (let yi = q - 1; yi >= 0; yi--) {
            for (let x = 0; x < q; x++) {
                // Find b such that y = a*x + b
                let b_val = 0;
                for (let b_try = 0; b_try < q; b_try++) {
                    if (field.add(field.mul(currentA, x), b_try) === yi) {
                        b_val = b_try;
                        break;
                    }
                }
                latinSquareEl.appendChild(makeCell(field.label(b_val), PALETTE[b_val % PALETTE.length], cellSize, cellFont));
            }
        }
    }
}

function makeText(content, x, y, cls, fontSize) {
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', x);
    t.setAttribute('y', y);
    t.setAttribute('class', cls);
    t.setAttribute('font-size', fontSize);
    t.textContent = content;
    return t;
}

function makeCell(text, color, size, fontSize) {
    const div = document.createElement('div');
    div.className = 'latin-cell';
    div.textContent = text;
    div.style.color = color;
    div.style.borderColor = color + '33';
    div.style.width = size + 'px';
    div.style.height = size + 'px';
    div.style.fontSize = fontSize + 'px';
    return div;
}

async function boot() {
    try {
        THREE = await import('three');
        const addons = await import('three/addons/controls/OrbitControls.js');
        OrbitControls = addons.OrbitControls;
        
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', init);
        } else {
            init();
        }
    } catch (e) {
        console.error("載入 Three.js 發生錯誤: ", e);
        alert("無法載入 3D 模組，請檢查網路連線或使用伺服器環境。");
    }
}

boot();
