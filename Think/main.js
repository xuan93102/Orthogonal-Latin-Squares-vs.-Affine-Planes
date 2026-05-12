let THREE, OrbitControls;

function formatAlphaLabel(coeffs) {
    const superscripts = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
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

function formatTupleLabel(coeffs) {
    return [...coeffs].reverse().join('');
}

function createGFPrime(p) {
    const labelAlpha = (i) => i.toString();
    const labelTuple = (i) => i.toString().padStart(2, '0');
    return {
        q: p,
        p,
        k: 1,
        add: (a, b) => (a + b) % p,
        mul: (a, b) => (a * b) % p,
        labelAlpha,
        labelTuple,
        label: labelAlpha,
        type: 'prime',
        desc: `質數體 GF(${p})\n加法: (a+b) mod ${p}\n乘法: (a×b) mod ${p}`
    };
}

function createGF2k(q, k) {
    // 不可約多項式 (含最高位):
    // k=1: x+1 (11_2 = 3) for GF(2)
    // k=2: x²+x+1 (111_2 = 7)
    // k=3: x³+x+1 (1011_2 = 11)
    // k=4: x⁴+x+1 (10011_2 = 19)
    // k=5: x⁵+x²+1 (100101_2 = 37)
    const polys = { 1: 3, 2: 7, 3: 11, 4: 19, 5: 37 };
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
        k,
        add: (a, b) => a ^ b,
        mul: mul,
        labelAlpha: (i) => formatAlphaLabel(intToBinaryCoeffs(i)),
        labelTuple: (i) => i.toString(2).padStart(k, '0'),
        label: (i) => formatAlphaLabel(intToBinaryCoeffs(i)),
        type: 'binary',
        desc: `擴張體 GF(2^${k}) = GF(${q})\n加法: 多項式 XOR (逐位元 XOR)\n乘法: 模不可約多項式算術`
    };
}

function createGFpk(p, k) {
    if (k === 1) return createGFPrime(p);
    if (p === 2) return createGF2k(Math.pow(2, k), k);

    const irreduciblePolys = {
        '3,2': [1, 0, 1], // x² + 1 over GF(3)
        '5,2': [2, 1, 1], // 2x² + x + 1 over GF(5)
        '7,2': [3, 1, 1], // 3x² + x + 1 over GF(7)
        '3,3': [1, 2, 0, 1], // x³ + 2x + 1 over GF(3)
        '5,3': [1, 0, 2, 1], // x³ + 2x + 1 over GF(5)
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

    const labelAlpha = (value) => formatAlphaLabel(intToCoeffs(value));
    const labelTuple = (value) => formatTupleLabel(intToCoeffs(value));

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
        k,
        add,
        mul,
        labelAlpha,
        labelTuple,
        label: labelAlpha,
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
    if (q === 32) return createGF2k(32, 5);

    const primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];
    if (primes.includes(q)) return createGFPrime(q);

    return null;
}

/** 加法反元素（小 q 暴力即可） */
function fieldNeg(field, x) {
    for (let i = 0; i < field.q; i++) {
        if (field.add(x, i) === 0) return i;
    }
    return 0;
}

/** 乘法反元素；0 則為 null */
function fieldInv(field, x) {
    if (x === 0) return null;
    for (let y = 1; y < field.q; y++) {
        if (field.mul(x, y) === 1) return y;
    }
    return null;
}

/** 2D：滿足 coefA·x + coefB·y = k 的所有格點（有序） */
function collectAffineLinePoints2D(field, coefA, coefB, k) {
    const q = field.q;
    const pts = [];
    if (coefA === 0 && coefB === 0) return pts;
    if (coefB !== 0) {
        const invB = fieldInv(field, coefB);
        for (let x = 0; x < q; x++) {
            const ax = field.mul(coefA, x);
            const rhs = field.add(k, fieldNeg(field, ax));
            const y = field.mul(invB, rhs);
            pts.push({ x, y });
        }
    } else {
        const invA = fieldInv(field, coefA);
        const x = field.mul(invA, k);
        for (let y = 0; y < q; y++) pts.push({ x, y });
    }
    pts.sort((p1, p2) => (p1.x !== p2.x ? p1.x - p2.x : p1.y - p2.y));
    return pts;
}

/** λ = ax + by + cz（格點上的線性式數值） */
function affineLambda(field, a, b, c, x, y, z) {
    return field.add(field.add(field.mul(a, x), field.mul(b, y)), field.mul(c, z));
}

/** 3D 高度：滿足 coefA·x + coefB·y + coefC·z = planeK 時，取 z 為垂直軸座標；coefC=0 時退回 z = ax+by+planeK */
function zFromAffinePlane(field, coefA, coefB, coefC, planeK, x, y) {
    if (coefC !== 0) {
        const invC = fieldInv(field, coefC);
        const ax = field.mul(coefA, x);
        const by = field.mul(coefB, y);
        const s = field.add(ax, by);
        const rhs = field.add(planeK, fieldNeg(field, s));
        return field.mul(invC, rhs);
    }
    const ax = field.mul(coefA, x);
    const by = field.mul(coefB, y);
    return field.add(field.add(ax, by), planeK);
}

const VALID_Q = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32];

const PALETTE = [
    '#5fa2ce', '#b48ead', '#d08770', '#a3be8c',
    '#ebcb8b', '#bf616a', '#88c0d0', '#81a1c1',
    '#f4a261', '#e76f51', '#2a9d8f', '#e9c46a',
    '#264653', '#6a4c93', '#c77dff', '#48cae4',
    '#f77f00', '#d62828', '#023e8a', '#80b918',
    '#bc6c25', '#606c38', '#ddb892', '#9b5de5',
    '#00bbf9',
];

// High-contrast colors for distinct plane identification (YZ, XZ, XY directions)
const PLANE_COLORS = [
    '#FF2D55', '#007AFF', '#34C759', '#FF9500',
    '#AF52DE', '#00C7BE', '#FF3B30', '#5856D6',
    '#FFCC00', '#FF6482', '#30D158', '#64D2FF',
    '#FFD60A', '#BF5AF2', '#FF453A', '#32D74B',
    '#0A84FF', '#FF6961', '#77DD77', '#FDFD96',
    '#AEC6CF', '#FFB347', '#B39EB5', '#CB99C9',
];

let currentField = null;
let currentA = 1;
let currentB = 1;
let currentC = 1;
let activeView = '2d';
let active2dPanel = 'line';
let active3dPanel = 'lines';
let showNumber = false;

let svgEl, latinSquareEl, eqDisplay2D, eqDisplay3D, aSlider, aSliderVal, fieldDesc;
let bSlider, bSliderVal;
let cSlider, cSliderVal;
let kInput, pkPInput, runKBtn;
let controls3d, showNumberCheck, view2d, view3d, view3dInteractive;
let interactiveVisualizer;
let interactiveCurveStyle = 'bezier';
let interactiveDisplayMode = 'lines'; // 'lines' | 'planes'
let axisXSelect, axisYSelect, axisZSelect;
let selectedAxisX = 'all', selectedAxisY = 'all', selectedAxisZ = 'all';
let valueSelect, valueSelectLabel, valueSelectWrapper, valuePlaceholder2D, valuePlaceholder3D;
let labelModeSelect;
let selectedValue = 'all';
let latinTooltip = null;
let curveStyle = 'bezier';   // 'straight' | 'smoothstep' | 'bezier' | 'catmull' | 'step' | 'arc'
let originMode = 'asc';      // 'asc' (0→α+1) | 'desc' (α+1→0)
let labelMode = 'algebraic'; // 'algebraic' | 'tuple'

let scene3dLines, scene3dSurface, scene3dLatin; // 3D 場景管理器執行個體 (kept for legacy compat)

// ===== Browser-tab 2D panel system =====
const PANEL_2D_DEFS = {
    line: { label: '線圖', title: '2D 線圖 (Lines)' },
    latin: { label: 'Latin方陣', title: 'Latin 方陣' },
};
const MAX_2D_PANELS = 4;
let open2DTabs = [];
let tab2DIdCounter = 0;
function create2DTabId() { return ++tab2DIdCounter; }

// ===== Browser-tab 3D panel system =====
const PANEL_DEFS = {
    lines: { label: '線圖', title: '3D 線圖 (Curves)', mode: '3d_lines' },
    linesPlanes: { label: '線圖＋疊層', title: '3D 線圖＋疊層平面', mode: '3d_lines_planes' },
    surface: { label: '3D線圖填色版', title: '3D 線圖填色版 (Curves)', mode: '3d_surface' },
    plotly: { label: '填色面', title: '3D 填色面 (Surface)', mode: 'plotly' },
    latin: { label: '疊層平面', title: '3D 疊層平面 (Latin Square)', mode: '3d_latin' },
};
const MAX_PANELS = 4;
let openTabs = [];      // { id, type, visualizer, plotlyEl, showNum, panelEl }
let tabIdCounter = 0;
function createTabId() { return ++tabIdCounter; }
let syncCameras = false;   // 是否同步所有視角
let syncSource = null;     // 正在被操作的 visualizer（同步來源）


class Visualizer3D {
    constructor(container) {
        this.container = container;
        this.bookmarks = [];
        this.lineBookmarks = [];
        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();
        this.init();
    }

    init() {
        this.scene = new THREE.Scene();
        const w = this.container.clientWidth || 400;
        const h = this.container.clientHeight || w;
        this.camera = new THREE.PerspectiveCamera(45, w / h, 0.1, 1000);
        this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
        this.renderer.setSize(w, h);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.camera.position.set(12, 12, 12);
        this.controls.update();
        this.lastQ = null;

        // 同步相機：當此 visualizer 被操作時，廣播給其他人
        this.controls.addEventListener('change', () => {
            if (syncCameras && syncSource === this) {
                broadcastCamera(this);
            }
        });
        this.renderer.domElement.addEventListener('pointerdown', (e) => {
            syncSource = this;
            this._pointerDownPos = { x: e.clientX, y: e.clientY };
        });
        this.renderer.domElement.addEventListener('pointerup', (e) => {
            // Only treat as click if pointer didn't move much (not a drag)
            if (!this._pointerDownPos) return;
            const dx = e.clientX - this._pointerDownPos.x;
            const dy = e.clientY - this._pointerDownPos.y;
            if (Math.sqrt(dx * dx + dy * dy) > 5) return; // was a drag

            const rect = this.renderer.domElement.getBoundingClientRect();
            this.mouse.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
            this.mouse.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

            this.raycaster.setFromCamera(this.mouse, this.camera);

            // 面模式：點擊書籤 Tab
            const bookmarkMeshes = this.bookmarks.map(b => b.mesh);
            const intersects = this.raycaster.intersectObjects(bookmarkMeshes);
            if (intersects.length > 0) {
                const clickedMesh = intersects[0].object;
                const bookmark = this.bookmarks.find(b => b.mesh === clickedMesh);
                if (bookmark) { showPlaneModal(bookmark); return; }
            }

            // 線模式：點擊 dot 節點
            if (this.lineBookmarks.length > 0) {
                const dotMeshes = this.lineBookmarks.map(b => b.mesh);
                const dotIntersects = this.raycaster.intersectObjects(dotMeshes);
                if (dotIntersects.length > 0) {
                    const clickedMesh = dotIntersects[0].object;
                    const lb = this.lineBookmarks.find(b => b.mesh === clickedMesh);
                    if (lb) showLineModal(lb);
                }
            }
        });

        const light = new THREE.DirectionalLight(0xffffff, 3);
        light.position.set(5, 10, 7.5);
        this.scene.add(light);
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.5));

        this.objects = new THREE.Group();
        this.scene.add(this.objects);

        this.animate = this.animate.bind(this);
        this.animate();

        // ResizeObserver：使用容器實際高度（若有），否則保持正方形
        this._ro = new ResizeObserver(() => {
            const w = this.container.clientWidth;
            if (w === 0) return;
            const h = this.container.clientHeight || w;
            const aspect = w / h;
            this.camera.aspect = aspect;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(w, h);
        });
        this._ro.observe(this.container);
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
        this.lineBookmarks = [];
    }

    render(mode, field, showText, aVal, bVal, cVal, filters = { x: 'all', y: 'all', z: 'all', c: 'all' }, panelCurveStyle) {
        this.clear();
        const q = field.q;
        const spacing = 1.6;
        const _curveStyle = panelCurveStyle || curveStyle;
        // coord: maps field index → world position
        // asc: 0 → 0, q-1 → (q-1)*spacing
        // desc: 0 → (q-1)*spacing, q-1 → 0
        const maxPos = (q - 1) * spacing;
        const coord = (i) => originMode === 'desc' ? (maxPos - i * spacing) : (i * spacing);
        const offset = 0; // no centering — axes start at 0

        if (this.lastQ !== q || this._lastOriginMode !== originMode) {
            const center = maxPos / 2;
            const distance = Math.max(16, q * spacing * 1.9);
            const camDist = this.container.id === 'view-3d-interactive-inner' ? distance * 1.33 : distance;
            if (originMode === 'desc') {
                // Y 軸向下模式：正面高俯角視角 (圖二風格)
                // Origin (field-0) is at world (maxPos, maxPos, maxPos)
                this.camera.position.set(
                    center - camDist * 1.1,
                    center + camDist * 0.6,
                    center - camDist * 1.1
                );
            } else {
                // Y 軸向上模式：標準對角線俯視（原始設定）
                this.camera.position.set(
                    center + camDist,
                    center + camDist * 0.8,
                    center + camDist
                );
            }
            this.camera.updateProjectionMatrix();
            this.controls.target.set(center, center, center);
            this.controls.update();
            this.lastQ = q;
            this._lastOriginMode = originMode;
        }

        this.drawAxes(q, spacing, coord, field, mode);

        const xFilter = filters.x === 'all' ? null : Number(filters.x);
        const yFilter = filters.y === 'all' ? null : Number(filters.y);
        const zFilter = filters.z === 'all' ? null : Number(filters.z);
        const cFilter = filters.c === 'all' ? null : Number(filters.c);

        if (mode === '3d_latin') {
            const group = new THREE.Group();
            for (let c = 0; c < q; c++) {
                if (cFilter !== null && c !== cFilter) continue;

                const color = PALETTE[c % PALETTE.length];
                const mat = new THREE.MeshPhongMaterial({
                    color: new THREE.Color(color), transparent: true, opacity: 0.65
                });

                for (let x = 0; x < q; x++) {
                    if (xFilter !== null && x !== xFilter) continue;
                    for (let y = 0; y < q; y++) {
                        if (yFilter !== null && y !== yFilter) continue;
                        const z = zFromAffinePlane(field, aVal, bVal, cVal, c, x, y);
                        if (zFilter !== null && z !== zFilter) continue;

                        const boxGeom = new THREE.BoxGeometry(0.08, 0.9, 0.9);
                        const mesh = new THREE.Mesh(boxGeom, mat);
                        mesh.position.set(coord(x), coord(z), coord(y));
                        group.add(mesh);

                        if (showText) {
                            const sprite = this.makeTextSprite(field.label(c), "white", 44);
                            sprite.position.set(coord(x) + 0.1, coord(z) + 0.5, coord(y));
                            sprite.scale.set(0.75, 0.75, 1);
                            group.add(sprite);
                        }
                    }
                }
            }
            this.objects.add(group);
        } else if (mode === '3d_lines' || mode === '3d_lines_planes') {
            if (mode === '3d_lines_planes') {
                this._drawStackedYZPlanes(field, aVal, bVal, cVal, xFilter, yFilter, zFilter, cFilter, coord, q, spacing, showText);
                // The user requested NO lines in this mode, so we skip drawing lines and dots.
            } else {
                this.lineBookmarks = [];
                for (let c = 0; c < q; c++) {
                    if (cFilter !== null && c !== cFilter) continue;

                    const color = PALETTE[c % PALETTE.length];
                    const group = new THREE.Group();
                    const material = new THREE.LineBasicMaterial({
                        color: color, transparent: true, opacity: 0.95, linewidth: 3
                    });
                    const dotMat = new THREE.MeshBasicMaterial({ color: color });
                    const dotGeom = new THREE.SphereGeometry(0.22, 12, 12);

                    const drawnDots = new Set();
                    const drawDot = (fx, fz, fy, lambdaVal) => {
                        const id = `${fx},${fz},${fy}`;
                        if (!drawnDots.has(id)) {
                            drawnDots.add(id);
                            const dot = new THREE.Mesh(dotGeom, dotMat);
                            dot.position.set(coord(fx), coord(fz), coord(fy));
                            group.add(dot);
                            this.lineBookmarks.push({
                                mesh: dot,
                                lambdaVal,
                                fieldX: fx,
                                fieldY: fy,
                                fieldZ: fz,
                                field, a: aVal, b: bVal, c: cVal, q,
                            });
                        }
                    };

                    for (let y = 0; y < q; y++) {
                        if (yFilter !== null && y !== yFilter) continue;
                        const points = [];
                        for (let x = 0; x < q; x++) {
                            if (xFilter !== null && x !== xFilter) continue;
                            const z = zFromAffinePlane(field, aVal, bVal, cVal, c, x, y);
                            if (zFilter !== null && z !== zFilter) continue;
                            points.push(new THREE.Vector3(coord(x), coord(z), coord(y)));
                            drawDot(x, z, y, c);
                        }
                        if (points.length <= 1) continue;
                        if (aVal === 0) {
                            group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
                        } else {
                            group.add(new THREE.Line(buildThreeCurveGeometry(points, points.length, 'x', _curveStyle), material));
                        }
                    }

                    for (let x = 0; x < q; x++) {
                        if (xFilter !== null && x !== xFilter) continue;
                        const points = [];
                        for (let y = 0; y < q; y++) {
                            if (yFilter !== null && y !== yFilter) continue;
                            const z = zFromAffinePlane(field, aVal, bVal, cVal, c, x, y);
                            if (zFilter !== null && z !== zFilter) continue;
                            points.push(new THREE.Vector3(coord(x), coord(z), coord(y)));
                            drawDot(x, z, y, c);
                        }
                        if (points.length <= 1) continue;
                        if (bVal === 0) {
                            group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
                        } else {
                            group.add(new THREE.Line(buildThreeCurveGeometry(points, points.length, 'y', _curveStyle), material));
                        }
                    }

                    for (let zf = 0; zf < q; zf++) {
                        if (zFilter !== null && zf !== zFilter) continue;
                        const rhs = field.add(c, fieldNeg(field, field.mul(cVal, zf)));
                        const points = [];
                        if (aVal === 0 && bVal === 0) continue;
                        if (bVal !== 0) {
                            const invB = fieldInv(field, bVal);
                            for (let x = 0; x < q; x++) {
                                if (xFilter !== null && x !== xFilter) continue;
                                const ax = field.mul(aVal, x);
                                const y = field.mul(invB, field.add(rhs, fieldNeg(field, ax)));
                                if (yFilter !== null && y !== yFilter) continue;
                                points.push(new THREE.Vector3(coord(x), coord(zf), coord(y)));
                                drawDot(x, zf, y, c);
                            }
                        } else {
                            const invA = fieldInv(field, aVal);
                            const x = field.mul(invA, rhs);
                            if (xFilter !== null && x !== xFilter) continue;
                            for (let y = 0; y < q; y++) {
                                if (yFilter !== null && y !== yFilter) continue;
                                points.push(new THREE.Vector3(coord(x), coord(zf), coord(y)));
                                drawDot(x, zf, y, c);
                            }
                        }
                        if (points.length <= 1) continue;
                        const isXConst = aVal !== 0 && bVal === 0;
                        if (isXConst) {
                            group.add(new THREE.Line(buildThreeCurveGeometry(points, points.length, 'y', _curveStyle), material));
                        } else if (aVal === 0) {
                            group.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
                        } else {
                            group.add(new THREE.Line(buildThreeCurveGeometry(points, points.length, 'x', _curveStyle), material));
                        }
                    }

                    this.objects.add(group);
                }
            }
        } else if (mode === '3d_surface') {
            const xAll = xFilter === null;
            const yAll = yFilter === null;
            const zAll = zFilter === null;
            const xRange = xAll ? Array.from({ length: q }, (_, i) => i) : [Number(xFilter)];
            const yRange = yAll ? Array.from({ length: q }, (_, i) => i) : [Number(yFilter)];

            const cRange = cFilter === null ? Array.from({ length: q }, (_, i) => i) : [Number(cFilter)];
            if (!xAll && yAll && zAll) {
                const x = Number(xFilter);
                const positions = [];
                const indices = [];
                for (let cIndex = 0; cIndex < cRange.length; cIndex++) {
                    const planeK = cRange[cIndex];
                    for (let y = 0; y < q; y++) {
                        const z = zFromAffinePlane(field, aVal, bVal, cVal, planeK, x, y);
                        positions.push(coord(x), coord(z), coord(y));
                    }
                }
                for (let cIndex = 0; cIndex < cRange.length - 1; cIndex++) {
                    for (let y = 0; y < q - 1; y++) {
                        const i = cIndex * q + y;
                        indices.push(i, i + 1, i + q);
                        indices.push(i + 1, i + q + 1, i + q);
                    }
                }
                const geom = new THREE.BufferGeometry();
                geom.setIndex(indices);
                geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
                geom.computeVertexNormals();
                this.objects.add(new THREE.Mesh(geom, new THREE.MeshPhongMaterial({
                    color: new THREE.Color(PALETTE[0]),
                    transparent: true, opacity: 0.55, side: THREE.DoubleSide, shininess: 30, flatShading: true
                })));
            } else if (xAll && !yAll && zAll) {
                const y = Number(yFilter);
                const positions = [];
                const indices = [];
                for (let cIndex = 0; cIndex < cRange.length; cIndex++) {
                    const planeK = cRange[cIndex];
                    for (let x = 0; x < q; x++) {
                        const z = zFromAffinePlane(field, aVal, bVal, cVal, planeK, x, y);
                        positions.push(coord(x), coord(z), coord(y));
                    }
                }
                for (let cIndex = 0; cIndex < cRange.length - 1; cIndex++) {
                    for (let x = 0; x < q - 1; x++) {
                        const i = cIndex * q + x;
                        indices.push(i, i + 1, i + q);
                        indices.push(i + 1, i + q + 1, i + q);
                    }
                }
                const geom = new THREE.BufferGeometry();
                geom.setIndex(indices);
                geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
                geom.computeVertexNormals();
                this.objects.add(new THREE.Mesh(geom, new THREE.MeshPhongMaterial({
                    color: new THREE.Color(PALETTE[1]),
                    transparent: true, opacity: 0.55, side: THREE.DoubleSide, shininess: 30, flatShading: true
                })));
            } else if (xAll && yAll && !zAll) {
                const z = coord(Number(zFilter));
                const positions = [];
                const indices = [];
                for (let y = 0; y < q; y++) {
                    for (let x = 0; x < q; x++) {
                        positions.push(coord(x), z, coord(y));
                    }
                }
                for (let y = 0; y < q - 1; y++) {
                    for (let x = 0; x < q - 1; x++) {
                        const i = y * q + x;
                        indices.push(i, i + 1, i + q);
                        indices.push(i + 1, i + q + 1, i + q);
                    }
                }
                const geom = new THREE.BufferGeometry();
                geom.setIndex(indices);
                geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
                geom.computeVertexNormals();
                this.objects.add(new THREE.Mesh(geom, new THREE.MeshPhongMaterial({
                    color: new THREE.Color('#8888ff'),
                    transparent: true, opacity: 0.25, side: THREE.DoubleSide, shininess: 10, flatShading: true
                })));
            } else {
                for (let c = 0; c < q; c++) {
                    if (cFilter !== null && c !== cFilter) continue;

                    const positions = [];
                    const indices = [];
                    for (let yIndex = 0; yIndex < yRange.length; yIndex++) {
                        const y = yRange[yIndex];
                        for (let xIndex = 0; xIndex < xRange.length; xIndex++) {
                            const x = xRange[xIndex];
                            const z = zFromAffinePlane(field, aVal, bVal, cVal, c, x, y);
                            if (!zAll && z !== Number(zFilter)) {
                                positions.push(NaN, NaN, NaN);
                            } else {
                                positions.push(coord(x), coord(z), coord(y));
                            }
                        }
                    }
                    if (xRange.length > 1 && yRange.length > 1) {
                        for (let y = 0; y < yRange.length - 1; y++) {
                            for (let x = 0; x < xRange.length - 1; x++) {
                                const i = y * xRange.length + x;
                                indices.push(i, i + 1, i + xRange.length);
                                indices.push(i + 1, i + xRange.length + 1, i + xRange.length);
                            }
                        }
                        const geom = new THREE.BufferGeometry();
                        geom.setIndex(indices);
                        geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
                        geom.computeVertexNormals();
                        this.objects.add(new THREE.Mesh(geom, new THREE.MeshPhongMaterial({
                            color: new THREE.Color(PALETTE[c % PALETTE.length]),
                            transparent: true,
                            opacity: 0.55,
                            side: THREE.DoubleSide,
                            shininess: 30,
                            flatShading: true
                        })));
                    } else if (xRange.length === 1 && yRange.length > 1) {
                        const positionsLine = [];
                        for (let i = 0; i < positions.length; i += 3) {
                            if (!Number.isNaN(positions[i])) {
                                positionsLine.push(new THREE.Vector3(positions[i], positions[i + 1], positions[i + 2]));
                            }
                        }
                        if (positionsLine.length > 0) {
                            const pointMat = new THREE.MeshPhongMaterial({ color: new THREE.Color(PALETTE[c % PALETTE.length]) });
                            const dotGeom = new THREE.SphereGeometry(0.12, 12, 12);
                            positionsLine.forEach(pt => {
                                const dot = new THREE.Mesh(dotGeom, pointMat);
                                dot.position.copy(pt);
                                this.objects.add(dot);
                            });
                        }
                    } else if (yRange.length === 1 && xRange.length > 1) {
                        const positionsLine = [];
                        for (let i = 0; i < positions.length; i += 3) {
                            if (!Number.isNaN(positions[i])) {
                                positionsLine.push(new THREE.Vector3(positions[i], positions[i + 1], positions[i + 2]));
                            }
                        }
                        if (positionsLine.length > 0) {
                            const pointMat = new THREE.MeshPhongMaterial({ color: new THREE.Color(PALETTE[c % PALETTE.length]) });
                            const dotGeom = new THREE.SphereGeometry(0.12, 12, 12);
                            positionsLine.forEach(pt => {
                                const dot = new THREE.Mesh(dotGeom, pointMat);
                                dot.position.copy(pt);
                                this.objects.add(dot);
                            });
                        }
                    }
                }
            }
        }
    }

    /** 沿 X、Y、Z 堆疊：大平面 + 凸出書籤 Tab + 格點著色 */
    _drawStackedYZPlanes(field, aVal, bVal, cVal, xFilter, yFilter, zFilter, cFilter, coord, q, spacing, showText) {
        const maxPos = (q - 1) * spacing;
        const group = new THREE.Group();
        this.bookmarks = [];

        const halfSp = spacing / 2;
        const tabW = spacing * 0.55;
        const tabH = spacing * 0.35;
        const tabDepth = 0.06;
        const tabGap = 0.02;
        const outerOffset = halfSp + tabH / 2 + tabGap;
        const xTabOffset = originMode === 'desc' ? maxPos + outerOffset : -outerOffset;
        const yTabOffset = originMode === 'desc' ? maxPos + outerOffset : -outerOffset;
        const zTabOffset = originMode === 'desc' ? maxPos + outerOffset : -outerOffset;

        const makeRoundedTabShape = (w, h, radius) => {
            const shape = new THREE.Shape();
            const r = Math.min(radius, w / 2, h / 2);
            shape.moveTo(-w / 2 + r, -h / 2);
            shape.lineTo(w / 2 - r, -h / 2);
            shape.quadraticCurveTo(w / 2, -h / 2, w / 2, -h / 2 + r);
            shape.lineTo(w / 2, h / 2 - r);
            shape.quadraticCurveTo(w / 2, h / 2, w / 2 - r, h / 2);
            shape.lineTo(-w / 2 + r, h / 2);
            shape.quadraticCurveTo(-w / 2, h / 2, -w / 2, h / 2 - r);
            shape.lineTo(-w / 2, -h / 2 + r);
            shape.quadraticCurveTo(-w / 2, -h / 2, -w / 2 + r, -h / 2);
            return shape;
        };

        const makeTab = (planeColor, pos, rotation, label, bookmarkData) => {
            const shape = makeRoundedTabShape(tabW, tabH, tabH * 0.35);
            const extrudeSettings = { depth: tabDepth, bevelEnabled: false };
            const tabGeom = new THREE.ExtrudeGeometry(shape, extrudeSettings);
            tabGeom.translate(0, 0, -tabDepth / 2);
            const tabMat = new THREE.MeshPhongMaterial({
                color: new THREE.Color(planeColor),
                transparent: true, opacity: 0.88, side: THREE.DoubleSide,
                emissive: new THREE.Color(planeColor), emissiveIntensity: 0.2,
            });
            const tabMesh = new THREE.Mesh(tabGeom, tabMat);
            tabMesh.position.copy(pos);
            if (rotation) tabMesh.rotation.copy(rotation);
            tabMesh.renderOrder = 100;
            group.add(tabMesh);

            const sprite = this.makeTextSprite(label, '#ffffff', 28);
            sprite.position.copy(pos);
            sprite.scale.set(0.55, 0.28, 1);
            sprite.renderOrder = 200;
            group.add(sprite);

            this.bookmarks.push({ mesh: tabMesh, ...bookmarkData });
        };

        const activeFilters = (xFilter !== null ? 1 : 0) + (yFilter !== null ? 1 : 0) + (zFilter !== null ? 1 : 0);
        // Show a plane type whenever its axis is filtered (so bookmarks appear for each filtered axis).
        // The isFullPlane check below prevents the opaque background from rendering when other filters are also active.
        const showYZ = activeFilters === 0 || xFilter !== null;
        const showXZ = activeFilters === 0 || yFilter !== null;
        const showXY = activeFilters === 0 || zFilter !== null;

        // ========== YZ 平面（固定 X）==========
        if (showYZ) {
            for (let x0 = 0; x0 < q; x0++) {
                if (xFilter !== null && x0 !== xFilter) continue;
                const ci = x0;
                const planeColor = PLANE_COLORS[ci % PLANE_COLORS.length];
                const isFullPlane = yFilter === null && zFilter === null;

                if (isFullPlane) {
                    const bgGeom = new THREE.PlaneGeometry(maxPos + spacing, maxPos + spacing);
                    const bgMat = new THREE.MeshPhongMaterial({
                        color: new THREE.Color(planeColor), transparent: false, opacity: 1.0,
                        side: THREE.DoubleSide, depthWrite: true,
                        polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1
                    });
                    const bgMesh = new THREE.Mesh(bgGeom, bgMat);
                    bgMesh.position.set(coord(x0), maxPos / 2, maxPos / 2);
                    bgMesh.rotation.y = Math.PI / 2;
                    bgMesh.renderOrder = -3;
                    group.add(bgMesh);
                }

                for (let y = 0; y < q; y++) {
                    if (yFilter !== null && y !== yFilter) continue;
                    for (let z = 0; z < q; z++) {
                        if (zFilter !== null && z !== zFilter) continue;
                        const lam = affineLambda(field, aVal, bVal, cVal, x0, y, z);
                        if (cFilter !== null && lam !== cFilter) continue;
                        const mat = new THREE.MeshPhongMaterial({
                            color: new THREE.Color(PALETTE[lam % PALETTE.length]),
                            transparent: true, opacity: 0.4, side: THREE.DoubleSide, depthWrite: false,
                        });
                        const geom = new THREE.PlaneGeometry(spacing * 0.9, spacing * 0.9);
                        const mesh = new THREE.Mesh(geom, mat);
                        mesh.position.set(coord(x0), coord(z), coord(y));
                        mesh.rotation.y = Math.PI / 2;
                        mesh.renderOrder = -2;
                        group.add(mesh);
                    }
                }

                const tabPos = new THREE.Vector3(coord(x0), yTabOffset, coord(0));
                const tabRot = new THREE.Euler(0, 0, 0);
                makeTab(planeColor, tabPos, tabRot, `X=${field.label(x0)}`, {
                    type: 'YZ', x: x0, y: null, z: null, field, a: aVal, b: bVal, c: cVal, q,
                });
            }
        }

        // ========== XZ 平面（固定 Y）==========
        if (showXZ) {
            for (let y0 = 0; y0 < q; y0++) {
                if (yFilter !== null && y0 !== yFilter) continue;
                const ci = q + y0;
                const planeColor = PLANE_COLORS[ci % PLANE_COLORS.length];
                const isFullPlane = xFilter === null && zFilter === null;

                if (isFullPlane) {
                    const bgGeom = new THREE.PlaneGeometry(maxPos + spacing, maxPos + spacing);
                    const bgMat = new THREE.MeshPhongMaterial({
                        color: new THREE.Color(planeColor), transparent: false, opacity: 1.0,
                        side: THREE.DoubleSide, depthWrite: true,
                        polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1
                    });
                    const bgMesh = new THREE.Mesh(bgGeom, bgMat);
                    bgMesh.position.set(maxPos / 2, maxPos / 2, coord(y0));
                    bgMesh.rotation.set(0, 0, 0);
                    bgMesh.renderOrder = -3;
                    group.add(bgMesh);
                }

                for (let x = 0; x < q; x++) {
                    if (xFilter !== null && x !== xFilter) continue;
                    for (let z = 0; z < q; z++) {
                        if (zFilter !== null && z !== zFilter) continue;
                        const lam = affineLambda(field, aVal, bVal, cVal, x, y0, z);
                        if (cFilter !== null && lam !== cFilter) continue;
                        const mat = new THREE.MeshPhongMaterial({
                            color: new THREE.Color(PALETTE[lam % PALETTE.length]),
                            transparent: true, opacity: 0.4, side: THREE.DoubleSide, depthWrite: false,
                        });
                        const geom = new THREE.PlaneGeometry(spacing * 0.9, spacing * 0.9);
                        const mesh = new THREE.Mesh(geom, mat);
                        mesh.position.set(coord(x), coord(z), coord(y0));
                        mesh.rotation.set(0, 0, 0);
                        mesh.renderOrder = -2;
                        group.add(mesh);
                    }
                }

                const tabPos = new THREE.Vector3(xTabOffset, coord(0), coord(y0));
                const tabRot = new THREE.Euler(originMode === 'desc' ? -Math.PI / 2 : Math.PI / 2, 0, 0);
                makeTab(planeColor, tabPos, tabRot, `Y=${field.label(y0)}`, {
                    type: 'XZ', x: null, y: y0, z: null, field, a: aVal, b: bVal, c: cVal, q,
                });
            }
        }

        // ========== XY 平面（固定 Z）==========
        if (showXY) {
            for (let z0 = 0; z0 < q; z0++) {
                if (zFilter !== null && z0 !== zFilter) continue;
                const ci = 2 * q + z0;
                const planeColor = PLANE_COLORS[ci % PLANE_COLORS.length];
                const isFullPlane = xFilter === null && yFilter === null;

                if (isFullPlane) {
                    const bgGeom = new THREE.PlaneGeometry(maxPos + spacing, maxPos + spacing);
                    const bgMat = new THREE.MeshPhongMaterial({
                        color: new THREE.Color(planeColor), transparent: false, opacity: 1.0,
                        side: THREE.DoubleSide, depthWrite: true,
                        polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1
                    });
                    const bgMesh = new THREE.Mesh(bgGeom, bgMat);
                    bgMesh.position.set(maxPos / 2, coord(z0), maxPos / 2);
                    bgMesh.rotation.x = Math.PI / 2;
                    bgMesh.renderOrder = -3;
                    group.add(bgMesh);
                }

                for (let x = 0; x < q; x++) {
                    if (xFilter !== null && x !== xFilter) continue;
                    for (let y = 0; y < q; y++) {
                        if (yFilter !== null && y !== yFilter) continue;
                        const lam = affineLambda(field, aVal, bVal, cVal, x, y, z0);
                        if (cFilter !== null && lam !== cFilter) continue;
                        const mat = new THREE.MeshPhongMaterial({
                            color: new THREE.Color(PALETTE[lam % PALETTE.length]),
                            transparent: true, opacity: 0.4, side: THREE.DoubleSide, depthWrite: false,
                        });
                        const geom = new THREE.PlaneGeometry(spacing * 0.9, spacing * 0.9);
                        const mesh = new THREE.Mesh(geom, mat);
                        mesh.position.set(coord(x), coord(z0), coord(y));
                        mesh.rotation.x = Math.PI / 2;
                        mesh.renderOrder = -2;
                        group.add(mesh);
                    }
                }

                const tabPos = new THREE.Vector3(coord(0), coord(z0), zTabOffset);
                const tabRot = new THREE.Euler(Math.PI / 2, 0, 0);
                makeTab(planeColor, tabPos, tabRot, `Z=${field.label(z0)}`, {
                    type: 'XY', x: null, y: null, z: z0, field, a: aVal, b: bVal, c: cVal, q,
                });
            }
        }

        this.objects.add(group);
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

    drawAxes(q, spacing, coord, field, mode) {
        const maxPos = (q - 1) * spacing;
        const isDesc = originMode === 'desc';
        const margin = 1.8;

        // The "origin corner" is where field-element 0 sits on ALL three axes.
        //   asc:  coord(0) = 0       → origin at world (0, 0, 0)  — left/bottom/back
        //   desc: coord(0) = maxPos  → origin at world (maxPos, maxPos, maxPos) — right/top/front
        //         BUT for the diagram we want origin at LEFT/TOP/BACK corner,
        //         meaning X and Z also start at coord(0)=maxPos and extend toward SMALLER world coords.
        const o = coord(0); // world coordinate of field-element 0 on every axis

        // Arrow directions:
        //   asc:  X→+x, Y→+y, Z→+z  (origin at small corner, extending outward)
        //   desc: X→-x, Y→-y, Z→-z  (origin at large corner, extending toward zero)
        const xDir = isDesc ? new THREE.Vector3(-1, 0, 0) : new THREE.Vector3(1, 0, 0);
        const yDir = isDesc ? new THREE.Vector3(0, -1, 0) : new THREE.Vector3(0, 1, 0);
        const zDir = isDesc ? new THREE.Vector3(0, 0, -1) : new THREE.Vector3(0, 0, 1);

        // Arrow start points: slightly behind the origin so the margin shows before tick 0
        const xStartOffset = isDesc ? margin : -margin;
        const zStartOffset = isDesc ? margin : -margin;

        const xStart = new THREE.Vector3(o + xStartOffset, o, o);
        const yStart = new THREE.Vector3(o, o, o);
        const zStart = new THREE.Vector3(o, o, o + zStartOffset);

        const arrowLen = maxPos + margin + 3;

        // X 軸 (紅色)
        this.objects.add(new THREE.ArrowHelper(xDir, xStart, arrowLen, 0xff4444));
        // Y 軸 (綠色)
        this.objects.add(new THREE.ArrowHelper(yDir, yStart, maxPos + 3, 0x44ff44));
        // Z 軸 (藍色)
        this.objects.add(new THREE.ArrowHelper(zDir, zStart, arrowLen, 0x4444ff));

        // --- Axis Titles (at arrow tips) ---
        const xTipX = isDesc ? o - maxPos - 3.5 : o + maxPos + 3.5;
        const yTipY = isDesc ? o - maxPos - 3.5 : o + maxPos + 3.5;
        const zTipZ = isDesc ? o - maxPos - 3.5 : o + maxPos + 3.5;

        const lxTitle = this.makeTextSprite('X', "#ff8888", 60);
        lxTitle.position.set(xTipX, o, o);
        lxTitle.scale.set(2, 2, 1);
        this.objects.add(lxTitle);

        // Physical Y axis (Green) corresponds to Mathematical Z
        const lyTitle = this.makeTextSprite('Z', "#88ff88", 60);
        lyTitle.position.set(o, yTipY, o);
        lyTitle.scale.set(2, 2, 1);
        this.objects.add(lyTitle);

        // Physical Z axis (Blue) corresponds to Mathematical Y
        const lzTitle = this.makeTextSprite('Y', "#8888ff", 60);
        lzTitle.position.set(o, o, zTipZ);
        lzTitle.scale.set(2, 2, 1);
        this.objects.add(lzTitle);

        // --- Tick labels ---
        // Labels sit just outside the origin corner on each axis.
        // "Outside" means away from the data cube, i.e., opposite to the arrow direction.
        const labelOff = 1.1;
        // Perpendicular offset for X and Z labels (on the Y axis):
        //   asc:  labels below origin → y = o - labelOff
        //   desc: labels above origin → y = o + labelOff
        const yLabelBase = isDesc ? o + labelOff : o - labelOff;
        // Perpendicular offset for Y labels (on the X axis):
        //   asc:  labels left of origin → x = o - labelOff
        //   desc: labels right of origin → x = o + labelOff
        const xLabelBase = isDesc ? o + labelOff : o - labelOff;

        for (let i = 0; i < q; i++) {
            const lbl = field.label(i);
            const size = 1.1;
            const ci = coord(i); // world position of field-element i

            // X tick: along x axis, offset in Y
            const lx = this.makeTextSprite(lbl, "#ffaaaa", 40);
            lx.position.set(ci, yLabelBase, o);
            lx.scale.set(size, size, 1);
            this.objects.add(lx);

            // Y tick: along y axis, offset in X
            const ly = this.makeTextSprite(lbl, "#aaffaa", 40);
            ly.position.set(xLabelBase, ci, o);
            ly.scale.set(size, size, 1);
            this.objects.add(ly);

            // Z tick: along z axis, offset in Y
            const lz = this.makeTextSprite(lbl, "#aaaaff", 40);
            lz.position.set(o, yLabelBase, ci);
            lz.scale.set(size, size, 1);
            this.objects.add(lz);
        }
    }

    makeTextSprite(message, color = null, fontSize = 42) {
        if (color === null) {
            color = document.body.classList.contains('light-mode') ? "#1a1a1a" : "white";
        }
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        const padding = 14;
        ctx.font = `bold ${fontSize}px sans-serif`;
        let width = Math.ceil(ctx.measureText(message).width) + padding * 2;
        let height = Math.ceil(fontSize * 1.6) + padding;

        if (width < 128) width = 128;
        if (height < 64) height = 64;

        if (width > 320) {
            const scale = 320 / width;
            fontSize = Math.max(16, Math.floor(fontSize * scale));
            ctx.font = `bold ${fontSize}px sans-serif`;
            width = Math.ceil(ctx.measureText(message).width) + padding * 2;
            height = Math.ceil(fontSize * 1.6) + padding;
        }

        canvas.width = width;
        canvas.height = height;
        ctx.fillStyle = color;
        ctx.font = `bold ${fontSize}px sans-serif`;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(message, width / 2, height / 2);

        const texture = new THREE.CanvasTexture(canvas);
        texture.minFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        const spriteMaterial = new THREE.SpriteMaterial({ map: texture, depthTest: false, transparent: true });
        const sprite = new THREE.Sprite(spriteMaterial);
        sprite.renderOrder = 999;
        const baseScale = 0.75;
        sprite.scale.set((width / 128) * baseScale, (height / 128) * baseScale, 1);
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
    cSlider = document.getElementById('c-slider');
    cSliderVal = document.getElementById('c-slider-val');
    fieldDesc = document.getElementById('field-desc');

    kInput = document.getElementById('k-input');
    pkPInput = document.getElementById('pk-p-input');
    runKBtn = document.getElementById('run-k');

    controls3d = document.getElementById('controls-3d');
    showNumberCheck = document.getElementById('show-number-check');
    axisXSelect = document.getElementById('axis-x-select');
    axisYSelect = document.getElementById('axis-y-select');
    axisZSelect = document.getElementById('axis-z-select');
    valueSelect = document.getElementById('value-select');
    valueSelectLabel = document.getElementById('value-select-label');
    valueSelectWrapper = document.getElementById('value-select-wrapper');
    valuePlaceholder2D = document.getElementById('value-placeholder-2d');
    valuePlaceholder3D = document.getElementById('value-placeholder-3d');
    labelModeSelect = document.getElementById('label-mode-select');
    view2d = document.getElementById('view-2d');
    view3d = document.getElementById('view-3d');
    view3dInteractive = document.getElementById('view-3d-interactive');

    // Init browser-tab 3D system — open first tab (lines) by default
    initBrowserTabs();
    addBrowserTab('lines');

    // Init browser-tab 2D system
    initBrowserTabs2D();

    interactiveVisualizer = new Visualizer3D(document.getElementById('view-3d-interactive-inner'));
    const interactiveOverlay = view3dInteractive.querySelector('.panel-curve-style-overlay');
    if (interactiveOverlay) {
        interactiveOverlay.querySelectorAll('.panel-curve-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                interactiveCurveStyle = btn.dataset.style;
                interactiveOverlay.querySelectorAll('.panel-curve-btn').forEach(b => b.classList.toggle('active', b.dataset.style === interactiveCurveStyle));
                if (activeView === '3d-interactive') render();
            });
        });
    }

    // 左上角：線/面 模式切換
    const displayModeOverlay = document.getElementById('interactive-display-mode-overlay');
    const lineHint = document.createElement('div');
    lineHint.className = 'interactive-line-hint';
    lineHint.id = 'interactive-line-hint';
    lineHint.textContent = '點擊節點查看平面資訊';
    const interactivePanelEl = view3dInteractive.querySelector('.browser-panel-lines-planes');
    if (interactivePanelEl) interactivePanelEl.appendChild(lineHint);

    const updateInteractiveMode = () => {
        if (displayModeOverlay) {
            displayModeOverlay.querySelectorAll('.panel-display-mode-btn').forEach(b =>
                b.classList.toggle('active', b.dataset.mode === interactiveDisplayMode)
            );
        }
        if (interactiveOverlay) {
            interactiveOverlay.classList.toggle('hidden-by-mode', interactiveDisplayMode === 'planes');
        }
        if (lineHint) lineHint.style.display = interactiveDisplayMode === 'lines' ? 'block' : 'none';
    };

    if (displayModeOverlay) {
        displayModeOverlay.querySelectorAll('.panel-display-mode-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                interactiveDisplayMode = btn.dataset.mode;
                updateInteractiveMode();
                if (activeView === '3d-interactive') render();
            });
        });
    }
    updateInteractiveMode();

    // 同步視角開關
    const syncCheck = document.getElementById('sync-cameras-check');
    if (syncCheck) {
        syncCheck.addEventListener('change', e => {
            syncCameras = e.target.checked;
            if (syncCameras && openTabs.length > 0) {
                // 找第一個 Three.js tab 作為同步來源
                const srcTab = openTabs.find(t => t.visualizer);
                if (srcTab) broadcastCamera(srcTab.visualizer);
            }
        });
    }

    latinTooltip = document.createElement('div');
    latinTooltip.className = 'latin-tooltip';
    document.body.appendChild(latinTooltip);

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
            alert(`k 必須在 1 到 ${maxK} 之間，且 q = p^k 不可超過 32`);
            return;
        }
        switchField(q, p, k);
    };

    pkPInput.addEventListener('change', () => updateKLimit());
    updateKLimit();

    showNumberCheck.onchange = e => {
        showNumber = e.target.checked;
        render();
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

    valueSelect.onchange = e => {
        selectedValue = e.target.value;
        updateValueSelectLabel();
        render();
        // Resize any open Plotly tabs
        if (activeView === '3d') {
            openTabs.filter(t => t.type === 'plotly').forEach(t => {
                if (t.plotlyEl && window.Plotly && Plotly.Plots && Plotly.Plots.resize) {
                    setTimeout(() => Plotly.Plots.resize(t.plotlyEl), 100);
                }
            });
        }
    };

    if (labelModeSelect) {
        labelModeSelect.value = labelMode;
        labelModeSelect.addEventListener('change', e => {
            labelMode = e.target.value === 'tuple' ? 'tuple' : 'algebraic';
            if (!currentField) return;
            applyFieldLabelMode();
            updateAxisSelectors();
            updateValueSelector();
            updateSliderLabelA();
            updateSliderLabelB();
            updateSliderLabelC();
            render();
        });
    }

    // Theme toggle
    const themeToggle = document.getElementById('theme-toggle-checkbox');
    const themeIconOverlay = document.getElementById('theme-icon-overlay');
    let themeIconTimeout = null;

    const showThemeIcon = icon => {
        if (!themeIconOverlay) return;
        themeIconOverlay.innerHTML = `<span>${icon}</span>`;
        themeIconOverlay.classList.add('active');
        clearTimeout(themeIconTimeout);
        themeIconTimeout = setTimeout(() => {
            themeIconOverlay.classList.remove('active');
        }, 550);
    };

    themeToggle.onchange = e => {
        const isLight = e.target.checked;
        document.body.classList.toggle('light-mode', isLight);
        const text = isLight ? '亮色' : '暗色';
        document.querySelector('.toggle-text').textContent = text;
        showThemeIcon(isLight ? '☀' : '🌙');
        // Re-render to update colors
        render();
        if (activeView === '3d') {
            openTabs.filter(t => t.type === 'plotly').forEach(t => {
                if (t.plotlyEl && window.Plotly && Plotly.Plots && Plotly.Plots.resize) {
                    setTimeout(() => Plotly.Plots.resize(t.plotlyEl), 100);
                }
            });
        }
    };

    const tab2D = document.getElementById('tab-2d');
    const tab3D = document.getElementById('tab-3d');
    const tab3DInteractive = document.getElementById('tab-3d-interactive');
    [tab2D, tab3D, tab3DInteractive].forEach(button => {
        if (!button) return;
        button.addEventListener('click', () => {
            if (button.id === 'tab-2d') setActiveView('2d');
            else if (button.id === 'tab-3d') setActiveView('3d');
            else if (button.id === 'tab-3d-interactive') setActiveView('3d-interactive');
        });
    });

    const tab2DLine = document.getElementById('tab-2d-line');
    const tab2DLatin = document.getElementById('tab-2d-latin');
    if (tab2DLine) tab2DLine.addEventListener('click', () => setActive2DPanel('line'));
    if (tab2DLatin) tab2DLatin.addEventListener('click', () => setActive2DPanel('latin'));

    // (3D panel tab management handled by browser-tab system)

    aSlider.addEventListener('input', e => {
        setA(parseInt(e.target.value, 10));
    });

    bSlider.addEventListener('input', e => {
        setB(parseInt(e.target.value, 10));
    });

    if (cSlider) {
        cSlider.addEventListener('input', e => {
            setC(parseInt(e.target.value, 10));
        });
    }

    // Curve style — panel buttons on 2D line chart (global, affects 2D only)
    document.querySelectorAll('.curve-panel-btn[data-style]').forEach(btn => {
        btn.addEventListener('click', () => {
            curveStyle = btn.dataset.style;
            document.querySelectorAll('.curve-panel-btn[data-style]').forEach(b => b.classList.toggle('active', b.dataset.style === curveStyle));
            render2D();
        });
    });

    // 3D Y-axis direction buttons (asc = arrow up / desc = arrow down)
    document.querySelectorAll('.y-axis-dir-btn[data-origin]').forEach(btn => {
        btn.addEventListener('click', () => {
            originMode = btn.dataset.origin;
            document.querySelectorAll('.y-axis-dir-btn[data-origin]').forEach(b => b.classList.toggle('active', b === btn));
            updateOriginBadge();
            render();
        });
    });

    switchField(4);
    add2DTab('line');
    updatePanelVisibility();
    updateOriginBadge();

    // 绑定模态框关闭按钮
    document.getElementById('plane-modal-close').addEventListener('click', closePlaneModal);
    document.getElementById('plane-modal').addEventListener('click', (e) => {
        if (e.target.id === 'plane-modal') closePlaneModal();
    });
}

let currentModalCurveStyle = 'bezier';
let currentModalCFilter = 'all'; // λ filter for the modal

/** 显示平面详情模态框 */
function showPlaneModal(bookmark) {
    const modal = document.getElementById('plane-modal');
    const titleEl = document.getElementById('plane-modal-title');
    const equationEl = document.getElementById('plane-modal-equation');
    const latinEl = document.getElementById('plane-modal-latin-square');
    const curveEl = document.getElementById('plane-modal-curve-svg');
    const styleRow = document.getElementById('plane-modal-curve-style-row');

    const field = bookmark.field;
    const q = bookmark.q;
    const a = bookmark.a;
    const b = bookmark.b;
    const c = bookmark.c;

    // 根据平面类型设置标题和方程式
    let title, equation, coef1, coef2;
    if (bookmark.type === 'YZ') {
        title = `YZ 平面 (X = ${field.label(bookmark.x)})`;
        coef1 = b; coef2 = c;
        equation = `λ = ${field.label(b)}·y + ${field.label(c)}·z`;
    } else if (bookmark.type === 'XZ') {
        title = `XZ 平面 (Y = ${field.label(bookmark.y)})`;
        coef1 = a; coef2 = c;
        equation = `λ = ${field.label(a)}·x + ${field.label(c)}·z`;
    } else if (bookmark.type === 'XY') {
        title = `XY 平面 (Z = ${field.label(bookmark.z)})`;
        coef1 = a; coef2 = b;
        equation = `λ = ${field.label(a)}·x + ${field.label(b)}·y`;
    }

    titleEl.textContent = title;
    equationEl.textContent = `方程式: ${equation}`;

    // 清空之前的内容
    latinEl.innerHTML = '';
    curveEl.innerHTML = '';
    currentModalCFilter = 'all'; // reset λ filter when opening new modal

    // 若從線模式點擊傳入預設 λ，自動預選
    if (bookmark._defaultLambda !== undefined && bookmark._defaultLambda !== null) {
        currentModalCFilter = String(bookmark._defaultLambda);
    }

    // ── λ value selector row ──────────────────────────────────────────────────
    let lambdaRow = document.getElementById('plane-modal-lambda-row');
    if (!lambdaRow) {
        lambdaRow = document.createElement('div');
        lambdaRow.id = 'plane-modal-lambda-row';
        lambdaRow.style.cssText = 'display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:4px;';
        // Insert before the grid
        const grid = document.querySelector('.plane-modal-grid');
        if (grid) grid.parentNode.insertBefore(lambdaRow, grid);
    }
    lambdaRow.innerHTML = '';
    const lambdaLabel = document.createElement('span');
    lambdaLabel.style.cssText = 'color:var(--accent-color);font-weight:600;font-size:0.93em;white-space:nowrap;';
    lambdaLabel.textContent = '顯示 λ =';
    lambdaRow.appendChild(lambdaLabel);
    const lambdaSlider = document.createElement('input');
    lambdaSlider.type = 'range';
    lambdaSlider.min = -1;
    lambdaSlider.max = q - 1;
    lambdaSlider.value = currentModalCFilter !== 'all' ? parseInt(currentModalCFilter) : -1;
    lambdaSlider.step = 1;
    lambdaSlider.style.cssText = 'width:160px;accent-color:var(--primary-color);cursor:pointer;';
    lambdaRow.appendChild(lambdaSlider);
    const lambdaVal = document.createElement('span');
    lambdaVal.style.cssText = 'font-family:monospace;font-size:1em;color:var(--primary-color);min-width:60px;background:var(--surface-2);border-radius:6px;padding:3px 8px;';
    lambdaVal.textContent = currentModalCFilter !== 'all' ? field.label(parseInt(currentModalCFilter)) : 'All';
    lambdaRow.appendChild(lambdaVal);
    // ── end λ selector ────────────────────────────────────────────────────────

    // Capture the global axis filters at modal-open time so the Latin Square reflects them
    const globalAxisFilters = { x: selectedAxisX, y: selectedAxisY, z: selectedAxisZ };

    const redrawAll = () => {
        latinEl.innerHTML = '';
        curveEl.innerHTML = '';
        const cf = currentModalCFilter === 'all' ? null : parseInt(currentModalCFilter);
        renderSliceLatinSquare(latinEl, field, q, bookmark.type, bookmark.x, bookmark.y, bookmark.z, a, b, c, cf, globalAxisFilters);
        render2DCurveForPlane(curveEl, field, bookmark.type, bookmark.x ?? bookmark.y ?? bookmark.z, a, b, c, q, currentModalCurveStyle, cf, globalAxisFilters);
    };

    lambdaSlider.oninput = () => {
        const v = parseInt(lambdaSlider.value);
        if (v === -1) {
            currentModalCFilter = 'all';
            lambdaVal.textContent = 'All';
        } else {
            currentModalCFilter = String(v);
            lambdaVal.textContent = `${field.label(v)}`;
        }
        redrawAll();
    };

    const redrawCurve = () => {
        const cf = currentModalCFilter === 'all' ? null : parseInt(currentModalCFilter);
        curveEl.innerHTML = '';
        render2DCurveForPlane(curveEl, field, bookmark.type, bookmark.x ?? bookmark.y ?? bookmark.z, a, b, c, q, currentModalCurveStyle, cf, globalAxisFilters);
    };

    // Render style buttons
    if (styleRow) {
        styleRow.innerHTML = '';
        const styles = [
            { id: 'bezier', label: 'Bézier' },
            { id: 'smoothstep', label: 'Smoothstep' },
            { id: 'catmull', label: 'Catmull-Rom' },
            { id: 'arc', label: '弧線' },
            { id: 'step', label: '階梯' },
            { id: 'line', label: '直線' }
        ];
        styles.forEach(st => {
            const btn = document.createElement('button');
            btn.className = `panel-btn curve-btn ${currentModalCurveStyle === st.id ? 'active' : ''}`;
            btn.textContent = st.label;
            btn.onclick = () => {
                currentModalCurveStyle = st.id;
                Array.from(styleRow.children).forEach(c => c.classList.remove('active'));
                btn.classList.add('active');
                redrawCurve();
            };
            styleRow.appendChild(btn);
        });
    }

    // Initial render — use redrawAll so Latin Square and curve are always in sync with currentModalCFilter
    redrawAll();

    modal.style.display = 'flex';
}

/** 渲染切片专用的 Latin Square（正確的 q×q 方陣）
 * @param {null|object} globalAxisFilters  { x, y, z } 字串('all'|'0'..'q-1')，用來 dim 不符合條件的格子
 */
function renderSliceLatinSquare(targetEl, field, q, planeType, xVal, yVal, zVal, coefA, coefB, coefC, cFilter, globalAxisFilters) {
    targetEl.innerHTML = '';

    const maxMatrixSize = Math.min(380, Math.max(240, q * 60));
    const gap = Math.max(3, Math.min(6, Math.floor(maxMatrixSize / q / 10)));
    const cellSize = Math.max(32, Math.min(70, Math.floor((maxMatrixSize - (q - 1) * gap) / q)));
    const cellFont = Math.max(9, Math.min(14, Math.floor(cellSize * 0.38)));

    targetEl.style.display = 'grid';
    targetEl.style.gap = `${gap}px`;
    targetEl.style.gridTemplateColumns = `repeat(${q}, ${cellSize}px)`;
    targetEl.style.width = `${q * cellSize + (q - 1) * gap + 20}px`;
    targetEl.style.height = `${q * cellSize + (q - 1) * gap + 20}px`;
    targetEl.style.padding = '10px';
    targetEl.style.background = 'var(--surface-2)';
    targetEl.style.borderRadius = '10px';

    // Derive per-axis numeric filters from globalAxisFilters
    const gxf = globalAxisFilters && globalAxisFilters.x !== 'all' ? parseInt(globalAxisFilters.x) : null;
    const gyf = globalAxisFilters && globalAxisFilters.y !== 'all' ? parseInt(globalAxisFilters.y) : null;
    const gzf = globalAxisFilters && globalAxisFilters.z !== 'all' ? parseInt(globalAxisFilters.z) : null;

    // The modal shows the full q×q slice plane. The only case to dim everything is when
    // the plane's own fixed axis is filtered to a different value (e.g. XZ plane fixes Y=0
    // but global filter has Y=1 selected — this plane doesn't exist in that view).
    let planeAxisBlocked = false;
    if (planeType === 'YZ' && gxf !== null && gxf !== xVal) planeAxisBlocked = true;
    if (planeType === 'XZ' && gyf !== null && gyf !== yVal) planeAxisBlocked = true;
    if (planeType === 'XY' && gzf !== null && gzf !== zVal) planeAxisBlocked = true;

    for (let row = 0; row < q; row++) {
        const rowIdx = q - 1 - row; // flip so 0 is at bottom
        for (let col = 0; col < q; col++) {
            let lam, detailText;
            // free1 = col axis value, free2 = rowIdx axis value (in terms of the plane's free axes)
            let freeAxis1Val = col, freeAxis2Val = rowIdx;

            if (planeType === 'YZ') {
                // free axes: y=col, z=rowIdx
                lam = affineLambda(field, coefA, coefB, coefC, xVal, col, rowIdx);
                detailText = `y=${field.label(col)} | z=${field.label(rowIdx)} | \u03bb=${field.label(lam)}`;
                // Check free axes against global filters
                var freeAxisBlocked = (gyf !== null && gyf !== col) || (gzf !== null && gzf !== rowIdx);
            } else if (planeType === 'XZ') {
                // free axes: x=col, z=rowIdx
                lam = affineLambda(field, coefA, coefB, coefC, col, yVal, rowIdx);
                detailText = `x=${field.label(col)} | z=${field.label(rowIdx)} | \u03bb=${field.label(lam)}`;
                var freeAxisBlocked = (gxf !== null && gxf !== col) || (gzf !== null && gzf !== rowIdx);
            } else {
                // XY plane: free axes: x=col, y=rowIdx
                lam = affineLambda(field, coefA, coefB, coefC, col, rowIdx, zVal);
                detailText = `x=${field.label(col)} | y=${field.label(rowIdx)} | \u03bb=${field.label(lam)}`;
                var freeAxisBlocked = (gxf !== null && gxf !== col) || (gyf !== null && gyf !== rowIdx);
            }

            // Dim when: the plane itself is blocked, a free axis doesn't match global filter, or λ filter doesn't match
            const lamMatched = (cFilter === null || cFilter === undefined) || lam === cFilter;
            const matched = !planeAxisBlocked && !freeAxisBlocked && lamMatched;
            const color = matched ? PALETTE[lam % PALETTE.length] : '#444';
            const cell = makeCell(field.label(lam), color, cellSize, cellFont, matched ? detailText : '');
            if (!matched) { cell.style.opacity = '0.18'; cell.style.color = '#555'; }
            targetEl.appendChild(cell);
        }
    }
}

/** 关闭平面详情模态框 */
function closePlaneModal() {
    const modal = document.getElementById('plane-modal');
    modal.style.display = 'none';
}

/** 線模式：點擊 dot 時顯示「該點所在 XZ 平面（固定 Y）」的詳細資訊 */
function showLineModal(lb) {
    const bookmark = {
        type: 'XZ',
        x: null,
        y: lb.fieldY,
        z: null,
        field: lb.field,
        a: lb.a,
        b: lb.b,
        c: lb.c,
        q: lb.q,
        _defaultLambda: lb.lambdaVal,
    };
    showPlaneModal(bookmark);
}

/** 绘制特定平面的 2D 曲线（按 λ 值分組繪製曲線） */
function render2DCurveForPlane(svgEl, field, planeType, fixedValue, a, b, c, q, curveStyle, cFilter, globalAxisFilters = null) {
    const padding = 28;
    const svgSize = 340;
    const plotSize = svgSize - padding * 2;
    const spacing = plotSize / (q - 1 || 1);
    const gxf = globalAxisFilters && globalAxisFilters.x !== 'all' ? parseInt(globalAxisFilters.x) : null;
    const gyf = globalAxisFilters && globalAxisFilters.y !== 'all' ? parseInt(globalAxisFilters.y) : null;
    const gzf = globalAxisFilters && globalAxisFilters.z !== 'all' ? parseInt(globalAxisFilters.z) : null;

    svgEl.setAttribute('viewBox', `0 0 ${svgSize} ${svgSize}`);
    svgEl.innerHTML = '';

    const px = (i) => padding + i * spacing;
    const py = (i) => padding + plotSize - i * spacing;

    // Draw grid
    const gridGroup = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    gridGroup.setAttribute('opacity', '0.25');
    for (let i = 0; i < q; i++) {
        gridGroup.appendChild(createSVGLine(px(i), py(0), px(i), py(q - 1), '0.5', 'var(--muted)'));
        gridGroup.appendChild(createSVGLine(px(0), py(i), px(q - 1), py(i), '0.5', 'var(--muted)'));
    }
    svgEl.appendChild(gridGroup);

    // Draw axes
    svgEl.appendChild(createSVGLine(px(0), py(0), px(q - 1), py(0), '1.5', 'var(--text-color)'));
    svgEl.appendChild(createSVGLine(px(0), py(0), px(0), py(q - 1), '1.5', 'var(--text-color)'));

    // Determine if the fixed axis of this plane is itself filtered out (whole plane dimmed)
    let planeAxisBlocked = false;
    if (planeType === 'YZ' && gxf !== null && gxf !== fixedValue) planeAxisBlocked = true;
    if (planeType === 'XZ' && gyf !== null && gyf !== fixedValue) planeAxisBlocked = true;
    if (planeType === 'XY' && gzf !== null && gzf !== fixedValue) planeAxisBlocked = true;

    // Collect ALL points grouped by λ value (always full q×q — modal shows entire plane)
    // Each point also stores whether it matches global free-axis filters
    const groups = {};
    for (let i1 = 0; i1 < q; i1++) {
        for (let i2 = 0; i2 < q; i2++) {
            let lam, ptX, ptY, ptFreeAxisBlocked;
            if (planeType === 'YZ') {
                lam = affineLambda(field, a, b, c, fixedValue, i1, i2);
                ptX = i1; ptY = i2;
                ptFreeAxisBlocked = (gyf !== null && gyf !== i1) || (gzf !== null && gzf !== i2);
            } else if (planeType === 'XZ') {
                lam = affineLambda(field, a, b, c, i1, fixedValue, i2);
                ptX = i1; ptY = i2;
                ptFreeAxisBlocked = (gxf !== null && gxf !== i1) || (gzf !== null && gzf !== i2);
            } else {
                lam = affineLambda(field, a, b, c, i1, i2, fixedValue);
                ptX = i1; ptY = i2;
                ptFreeAxisBlocked = (gxf !== null && gxf !== i1) || (gyf !== null && gyf !== i2);
            }
            if (!groups[lam]) groups[lam] = [];
            groups[lam].push({ x: ptX, y: ptY, freeAxisBlocked: ptFreeAxisBlocked });
        }
    }

    // Draw curves for each λ group
    for (const [lamStr, pts] of Object.entries(groups)) {
        const lam = parseInt(lamStr);
        if (cFilter !== null && cFilter !== undefined && lam !== cFilter) continue;
        const color = PALETTE[lam % PALETTE.length];

        // Dim the whole group if the plane's fixed axis is blocked
        const groupDimmed = planeAxisBlocked;

        // Separate matched vs dimmed points for correct curve rendering
        const matchedPts = pts.filter(p => !p.freeAxisBlocked);
        const dimmedPts = pts.filter(p => p.freeAxisBlocked);

        // Helper to draw a path through a sorted point set
        const drawCurveSegment = (sortedPts, dimmed) => {
            if (sortedPts.length < 2) return;
            const isVert = sortedPts.every(p => p.x === sortedPts[0].x);
            let d = `M ${px(sortedPts[0].x)} ${py(sortedPts[0].y)}`;
            for (let k = 1; k < sortedPts.length; k++) {
                const x0 = px(sortedPts[k - 1].x), y0 = py(sortedPts[k - 1].y);
                const x1 = px(sortedPts[k].x), y1 = py(sortedPts[k].y);
                if (curveStyle === 'bezier') {
                    const dx = isVert ? 0 : (x1 - x0) * 0.42;
                    const dy = isVert ? (y1 - y0) * 0.42 : 0;
                    d += ` C ${x0 + dx} ${y0 + dy}, ${x1 - dx} ${y1 - dy}, ${x1} ${y1}`;
                } else if (curveStyle === 'smoothstep') {
                    const dx = isVert ? 0 : (x1 - x0) * 0.33;
                    const dy = isVert ? (y1 - y0) * 0.33 : 0;
                    d += ` C ${x0 + dx} ${y0 + dy}, ${x1 - dx} ${y1 - dy}, ${x1} ${y1}`;
                } else if (curveStyle === 'catmull') {
                    const p_prev = k - 2 >= 0 ? sortedPts[k - 2] : sortedPts[k - 1];
                    const p_next = k + 1 < sortedPts.length ? sortedPts[k + 1] : sortedPts[k];
                    const px_prev = px(p_prev.x), py_prev = py(p_prev.y);
                    const px_next = px(p_next.x), py_next = py(p_next.y);
                    const cp1x = x0 + (x1 - px_prev) / 6;
                    const cp1y = y0 + (y1 - py_prev) / 6;
                    const cp2x = x1 - (px_next - x0) / 6;
                    const cp2y = y1 - (py_next - y0) / 6;
                    d += ` C ${cp1x} ${cp1y}, ${cp2x} ${cp2y}, ${x1} ${y1}`;
                } else if (curveStyle === 'step') {
                    if (isVert || y0 === y1) { d += ` L ${x1} ${y1}`; }
                    else { const midX = x0 + (x1 - x0) / 2; d += ` L ${midX} ${y0} L ${midX} ${y1} L ${x1} ${y1}`; }
                } else if (curveStyle === 'arc') {
                    if (isVert) { d += ` L ${x1} ${y1}`; }
                    else { const midX = (x0 + x1) / 2; const midY = (y0 + y1) / 2 - Math.abs(x1 - x0) * 0.4; d += ` Q ${midX} ${midY}, ${x1} ${y1}`; }
                } else { d += ` L ${x1} ${y1}`; }
            }
            const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
            path.setAttribute('d', d);
            path.setAttribute('stroke', dimmed ? '#555' : color);
            path.setAttribute('stroke-width', dimmed ? '1.5' : '2.5');
            path.setAttribute('fill', 'none');
            path.setAttribute('opacity', dimmed ? '0.13' : '0.85');
            svgEl.appendChild(path);
        };

        if (groupDimmed) {
            // Whole group is dimmed (fixed axis blocked): draw all as one dimmed curve
            const isVert = pts.every(p => p.x === pts[0].x);
            const sortedAll = [...pts].sort((a, b) => isVert ? a.y - b.y : a.x - b.x);
            drawCurveSegment(sortedAll, true);
        } else {
            // Draw matched and dimmed segments separately
            const isVertM = matchedPts.length > 0 && matchedPts.every(p => p.x === matchedPts[0].x);
            const sortedMatched = [...matchedPts].sort((a, b) => isVertM ? a.y - b.y : a.x - b.x);
            drawCurveSegment(sortedMatched, false);

            const isVertD = dimmedPts.length > 0 && dimmedPts.every(p => p.x === dimmedPts[0].x);
            const sortedDimmed = [...dimmedPts].sort((a, b) => isVertD ? a.y - b.y : a.x - b.x);
            drawCurveSegment(sortedDimmed, true);
        }

        // Draw dots (per-point dim based on free-axis filter match)
        pts.forEach(p => {
            const ptDimmed = groupDimmed || p.freeAxisBlocked;
            const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
            circle.setAttribute('cx', px(p.x));
            circle.setAttribute('cy', py(p.y));
            circle.setAttribute('r', ptDimmed ? '3.5' : '4.5');
            circle.setAttribute('fill', ptDimmed ? '#555' : color);
            circle.setAttribute('stroke', 'var(--surface-1)');
            circle.setAttribute('stroke-width', '1.5');
            circle.setAttribute('opacity', ptDimmed ? '0.13' : '1');
            svgEl.appendChild(circle);
        });
    }

    // Axis labels
    for (let i = 0; i < q; i++) {
        const lbl = field.label(i);
        const txBottom = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        txBottom.setAttribute('x', px(i));
        txBottom.setAttribute('y', svgSize - 5);
        txBottom.setAttribute('text-anchor', 'middle');
        txBottom.setAttribute('fill', 'var(--muted)');
        txBottom.setAttribute('font-size', '10');
        txBottom.textContent = lbl;
        svgEl.appendChild(txBottom);

        const txLeft = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        txLeft.setAttribute('x', padding - 8);
        txLeft.setAttribute('y', py(i) + 3);
        txLeft.setAttribute('text-anchor', 'end');
        txLeft.setAttribute('fill', 'var(--muted)');
        txLeft.setAttribute('font-size', '10');
        txLeft.textContent = lbl;
        svgEl.appendChild(txLeft);
    }
}

function createSVGLine(x1, y1, x2, y2, strokeWidth, stroke) {
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', x1);
    line.setAttribute('y1', y1);
    line.setAttribute('x2', x2);
    line.setAttribute('y2', y2);
    line.setAttribute('stroke-width', strokeWidth);
    line.setAttribute('stroke', stroke);
    return line;
}

function updateOriginBadge() {
    const badge = document.getElementById('origin-mode-badge');
    if (!badge) return;
    const labels = { asc: 'Y軸: 0 在下方 ↑', desc: 'Y軸: 0 在上方 ↓' };
    badge.textContent = labels[originMode] || '';
}

// Build SVG path data for a set of points using the selected curveStyle
function buildSVGPath(pts, px, py, style) {
    style = style || curveStyle;
    if (pts.length === 0) return '';
    let d = `M ${px(pts[0].x)} ${py(pts[0].y)}`;
    if (style === 'straight' || pts.length === 1) {
        for (let k = 1; k < pts.length; k++) d += ` L ${px(pts[k].x)} ${py(pts[k].y)}`;
    } else if (style === 'smoothstep') {
        for (let k = 0; k < pts.length - 1; k++) {
            const x0 = px(pts[k].x), y0 = py(pts[k].y);
            const x1 = px(pts[k + 1].x), y1 = py(pts[k + 1].y);
            const dx = x1 - x0;
            const cp1x = x0 + dx * 0.33, cp2x = x1 - dx * 0.33;
            d += ` C ${cp1x} ${y0}, ${cp2x} ${y1}, ${x1} ${y1}`;
        }
    } else if (style === 'catmull') {
        // Catmull-Rom spline: uses neighboring points as control tangents
        const screenPts = pts.map(p => ({ x: px(p.x), y: py(p.y) }));
        for (let k = 0; k < screenPts.length - 1; k++) {
            const p0 = screenPts[Math.max(0, k - 1)];
            const p1 = screenPts[k];
            const p2 = screenPts[k + 1];
            const p3 = screenPts[Math.min(screenPts.length - 1, k + 2)];
            const cp1x = p1.x + (p2.x - p0.x) / 6;
            const cp1y = p1.y + (p2.y - p0.y) / 6;
            const cp2x = p2.x - (p3.x - p1.x) / 6;
            const cp2y = p2.y - (p3.y - p1.y) / 6;
            d += ` C ${cp1x} ${cp1y}, ${cp2x} ${cp2y}, ${p2.x} ${p2.y}`;
        }
    } else if (style === 'step') {
        // Step function: horizontal then vertical
        for (let k = 0; k < pts.length - 1; k++) {
            const x0 = px(pts[k].x), y0 = py(pts[k].y);
            const x1 = px(pts[k + 1].x), y1 = py(pts[k + 1].y);
            d += ` H ${x1} V ${y1}`;
        }
    } else if (style === 'arc') {
        // Arc: SVG arc segments between consecutive points
        for (let k = 0; k < pts.length - 1; k++) {
            const x0 = px(pts[k].x), y0 = py(pts[k].y);
            const x1 = px(pts[k + 1].x), y1 = py(pts[k + 1].y);
            const dx = x1 - x0, dy = y1 - y0;
            const r = Math.sqrt(dx * dx + dy * dy) * 0.65;
            const sweep = (dy < 0) ? 1 : 0;
            d += ` A ${r} ${r} 0 0 ${sweep} ${x1} ${y1}`;
        }
    } else {
        // bezier (default): S-curve between each pair
        for (let k = 0; k < pts.length - 1; k++) {
            const x0 = px(pts[k].x), y0 = py(pts[k].y);
            const x1 = px(pts[k + 1].x), y1 = py(pts[k + 1].y);
            const cpd = (x1 - x0) * 0.42;
            d += ` C ${x0 + cpd} ${y0}, ${x1 - cpd} ${y1}, ${x1} ${y1}`;
        }
    }
    return d;
}

// Build Three.js BufferGeometry for a polyline using the given style
function buildThreeCurveGeometry(points, q, direction, style) {
    style = style || curveStyle;
    if (style === 'straight') {
        return new THREE.BufferGeometry().setFromPoints(points);
    } else if (style === 'smoothstep') {
        const curvePath = new THREE.CurvePath();
        for (let k = 0; k < points.length - 1; k++) {
            const p0 = points[k], p1 = points[k + 1];
            const dx = direction === 'x' ? (p1.x - p0.x) * 0.33 : 0;
            const dz = direction === 'y' ? (p1.z - p0.z) * 0.33 : 0;
            const cp1 = new THREE.Vector3(p0.x + dx, p0.y, p0.z + dz);
            const cp2 = new THREE.Vector3(p1.x - dx, p1.y, p1.z - dz);
            curvePath.add(new THREE.CubicBezierCurve3(p0, cp1, cp2, p1));
        }
        return new THREE.BufferGeometry().setFromPoints(curvePath.getPoints(q * 8));
    } else if (style === 'catmull') {
        // Catmull-Rom via CatmullRomCurve3
        const curve = new THREE.CatmullRomCurve3(points);
        return new THREE.BufferGeometry().setFromPoints(curve.getPoints(q * 8));
    } else if (style === 'step') {
        // Step function: horizontal then vertical in 3D
        const stepPts = [];
        for (let k = 0; k < points.length - 1; k++) {
            const p0 = points[k], p1 = points[k + 1];
            stepPts.push(p0);
            if (direction === 'x') {
                stepPts.push(new THREE.Vector3(p1.x, p0.y, p0.z));
            } else {
                stepPts.push(new THREE.Vector3(p0.x, p0.y, p1.z));
            }
            stepPts.push(p1);
        }
        return new THREE.BufferGeometry().setFromPoints(stepPts);
    } else if (style === 'arc') {
        // Arc: quadratic bezier with apex lifted in y
        const curvePath = new THREE.CurvePath();
        for (let k = 0; k < points.length - 1; k++) {
            const p0 = points[k], p1 = points[k + 1];
            const mid = new THREE.Vector3(
                (p0.x + p1.x) / 2,
                (p0.y + p1.y) / 2 + Math.abs(p1.y - p0.y) * 0.5 + 0.6,
                (p0.z + p1.z) / 2
            );
            curvePath.add(new THREE.QuadraticBezierCurve3(p0, mid, p1));
        }
        return new THREE.BufferGeometry().setFromPoints(curvePath.getPoints(q * 8));
    } else {
        // bezier: existing S-curve logic
        const curvePath = new THREE.CurvePath();
        for (let k = 0; k < points.length - 1; k++) {
            const p0 = points[k], p1 = points[k + 1];
            const dx = direction === 'x' ? (p1.x - p0.x) * 0.42 : 0;
            const dz = direction === 'y' ? (p1.z - p0.z) * 0.42 : 0;
            const cp1 = new THREE.Vector3(p0.x + dx, p0.y, p0.z + dz);
            const cp2 = new THREE.Vector3(p1.x - dx, p1.y, p1.z - dz);
            curvePath.add(new THREE.CubicBezierCurve3(p0, cp1, cp2, p1));
        }
        return new THREE.BufferGeometry().setFromPoints(curvePath.getPoints(q * 8));
    }
}

function isPrime(n) {
    if (n < 2) return false;
    for (let i = 2; i <= Math.sqrt(n); i++) if (n % i === 0) return false;
    return true;
}

function getMaxKForPrime(p) {
    return Math.max(1, Math.floor(Math.log(32) / Math.log(p)));
}

function updateKLimit() {
    if (!pkPInput) return;
    const p = parseInt(pkPInput.value);
    const maxK = getMaxKForPrime(p);
    kInput.max = maxK;
    if (parseInt(kInput.value) > maxK) kInput.value = maxK;
}

function setActiveView(mode) {
    activeView = mode;
    const tab2D = document.getElementById('tab-2d');
    const tab3D = document.getElementById('tab-3d');
    const tab3DInteractive = document.getElementById('tab-3d-interactive');
    if (tab2D) {
        tab2D.classList.toggle('active', mode === '2d');
        tab2D.setAttribute('aria-selected', mode === '2d');
    }
    if (tab3D) {
        tab3D.classList.toggle('active', mode === '3d');
        tab3D.setAttribute('aria-selected', mode === '3d');
    }
    if (tab3DInteractive) {
        tab3DInteractive.classList.toggle('active', mode === '3d-interactive');
        tab3DInteractive.setAttribute('aria-selected', mode === '3d-interactive');
    }
    updateViewMode();
}

function updateViewMode() {
    const is3d = activeView === '3d';
    const is3dInteractive = activeView === '3d-interactive';
    const isAny3d = is3d || is3dInteractive;

    if (view2d) view2d.style.display = activeView === '2d' ? 'flex' : 'none';
    if (view3d) view3d.style.display = is3d ? 'flex' : 'none';
    if (view3dInteractive) view3dInteractive.style.display = is3dInteractive ? 'flex' : 'none';
    if (controls3d) controls3d.style.display = isAny3d ? 'flex' : 'none';

    const cSliderGroup = document.getElementById('c-slider-group');
    if (cSliderGroup) cSliderGroup.style.display = isAny3d ? 'flex' : 'none';

    document.getElementById('a-slider-label-text').innerText = '選擇係數 a（x）';
    updateValueControlPlacement();
    aSlider.min = 0;
    if (eqDisplay3D) eqDisplay3D.style.display = isAny3d ? 'block' : 'none';
    if (eqDisplay2D) eqDisplay2D.style.display = !isAny3d ? 'block' : 'none';

    if (isAny3d) {
        setTimeout(() => window.dispatchEvent(new Event('resize')), 10);
    }
    updatePanelVisibility();
    render();
}

function setActive2DPanel(panel) {
    active2dPanel = panel;
    updatePanelVisibility();
}

function applyFieldLabelMode() {
    if (!currentField) return;
    currentField.label = labelMode === 'tuple' ? currentField.labelTuple : currentField.labelAlpha;
}

function updateValueControlPlacement() {
    if (!valueSelectWrapper || !valuePlaceholder2D || !valuePlaceholder3D) return;
    if (activeView === '3d' || activeView === '3d-interactive') {
        valuePlaceholder3D.appendChild(valueSelectWrapper);
    } else {
        valuePlaceholder2D.appendChild(valueSelectWrapper);
    }
}

function setActive3DPanel(panel) {
    // Legacy: no-op — panel management now done via browser-tab system
    active3dPanel = panel;
}

function updatePanelVisibility() {
    // 3D panels are handled by browser-tab system — no classic panel-content toggling needed
}

function switchField(q, p = null, k = null) {
    currentField = createField(q, p, k);
    if (!currentField) return;
    applyFieldLabelMode();

    currentA = 1;
    currentB = 1;
    currentC = 1;
    selectedAxisX = 'all';
    selectedAxisY = 'all';
    selectedAxisZ = 'all';
    selectedValue = 'all';

    aSlider.min = 0;
    aSlider.max = q - 1;
    aSlider.value = currentA;

    updateAxisSelectors();
    updateValueSelector();

    bSlider.min = 0;
    bSlider.max = q - 1;
    bSlider.value = currentB;

    if (cSlider) {
        cSlider.min = 0;
        cSlider.max = q - 1;
        cSlider.value = currentC;
    }

    fieldDesc.innerText = currentField.desc;
    updateSliderLabelA();
    updateSliderLabelB();
    updateSliderLabelC();
    updateViewMode();
}

function setA(val) {
    currentA = parseInt(val, 10);
    aSlider.value = currentA;
    updateSliderLabelA();
    render();
}

function setB(val) {
    currentB = parseInt(val, 10);
    bSlider.value = currentB;
    updateSliderLabelB();
    render();
}

function setC(val) {
    currentC = parseInt(val, 10);
    if (cSlider) cSlider.value = currentC;
    updateSliderLabelC();
    render();
}

function updateAxisSelectors() {
    if (!currentField || !axisXSelect || !axisYSelect || !axisZSelect) return;
    const q = currentField.q;

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
            const label = currentField.label(i);
            select.appendChild(makeOption(String(i), label));
        }
    });

    axisXSelect.value = selectedAxisX;
    axisYSelect.value = selectedAxisY;
    axisZSelect.value = selectedAxisZ;
}

function updateValueSelector() {
    if (!currentField || !valueSelect) return;
    const q = currentField.q;

    const makeOption = (value, label) => {
        const opt = document.createElement('option');
        opt.value = value;
        opt.textContent = label;
        return opt;
    };

    valueSelect.innerHTML = '';
    valueSelect.appendChild(makeOption('all', 'All'));
    for (let i = 0; i < q; i++) {
        const label = currentField.label(i);
        valueSelect.appendChild(makeOption(String(i), label));
    }

    valueSelect.value = selectedValue;
    updateValueSelectLabel();
}

function updateValueSelectLabel() {
    if (selectedValue === 'all') {
        valueSelectLabel.textContent = '顯示: 全部';
    } else {
        valueSelectLabel.textContent = `顯示: ${currentField.label(parseInt(selectedValue))}`;
    }
}

function updateSliderLabelA() {
    if (!currentField || !aSliderVal) return;
    aSliderVal.textContent = `a = ${currentField.label(currentA)}`;
}

function updateSliderLabelB() {
    if (!currentField || !bSliderVal) return;
    bSliderVal.textContent = `b = ${currentField.label(currentB)}`;
}

function updateSliderLabelC() {
    if (!currentField || !cSliderVal) return;
    cSliderVal.textContent = `c = ${currentField.label(currentC)}`;
}

function render() {
    if (activeView === '3d' || activeView === '3d-interactive') {
        eqDisplay3D.innerHTML = `3D: λ = ${currentField.label(currentA)}·x + ${currentField.label(currentB)}·y + ${currentField.label(currentC)}·z`;
        if (activeView === '3d') {
            renderAllOpenTabs();
        } else if (activeView === '3d-interactive' && interactiveVisualizer) {
            const axisFilters = {
                x: selectedAxisX,
                y: selectedAxisY,
                z: selectedAxisZ,
                c: selectedValue
            };
            const interactiveMode = interactiveDisplayMode === 'planes' ? '3d_lines_planes' : '3d_lines';
            interactiveVisualizer.render(interactiveMode, currentField, false, currentA, currentB, currentC, axisFilters, interactiveCurveStyle);
        }
    }
    render2D();
    renderLatinSquare(currentField, currentField.q);
    renderAllOpen2DTabs();
}

function animate3DRedraw() {
    // Handled per-tab via Visualizer3D.animate loop
}

function renderPlotlySurface(field, aVal, bVal, cVal, filters, g) {
    if (!window.Plotly) return;
    if (!g) {
        // legacy fallback
        g = document.getElementById('viz-3d-plotly');
    }
    if (!g) return;

    const isLightMode = document.body.classList.contains('light-mode');
    const axisColor = isLightMode ? '#1a73e8' : '#ff4444';
    const axisColorG = isLightMode ? '#34a853' : '#44ff44';
    const axisColorB = isLightMode ? '#ea4335' : '#4444ff';
    const textColor = isLightMode ? '#202124' : '#ffaaaa';
    const textColorG = isLightMode ? '#202124' : '#aaffaa';
    const textColorB = isLightMode ? '#202124' : '#aaaaff';

    const q = field.q;
    const cFilter = filters.c === 'all' ? null : parseInt(filters.c);

    const xVals = [];
    const yVals = [];
    const tickTexts = [];
    for (let i = 0; i < q; i++) {
        xVals.push(i);
        yVals.push(i);
        tickTexts.push(field.label(i));
    }

    const data = [];
    const colorscales = ['Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis', 'Blues', 'Reds', 'Greens'];
    const maxVal = q - 1;

    data.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: [0, maxVal],
        y: [0, 0],
        z: [0, 0],
        line: { color: axisColor, width: 6 },
        marker: { size: 3, color: axisColor },
        showlegend: false,
        hoverinfo: 'skip'
    });
    data.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: [0, 0],
        y: [0, maxVal],
        z: [0, 0],
        line: { color: axisColorG, width: 6 },
        marker: { size: 3, color: axisColorG },
        showlegend: false,
        hoverinfo: 'skip'
    });
    data.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: [0, 0],
        y: [0, 0],
        z: [0, maxVal],
        line: { color: axisColorB, width: 6 },
        marker: { size: 3, color: axisColorB },
        showlegend: false,
        hoverinfo: 'skip'
    });

    for (let c = 0; c < q; c++) {
        if (cFilter !== null && c !== cFilter) continue;

        const zData = [];
        for (let y = 0; y < q; y++) {
            const row = [];
            for (let x = 0; x < q; x++) {
                row.push(zFromAffinePlane(field, aVal, bVal, cVal, c, x, y));
            }
            zData.push(row);
        }

        data.push({
            z: zData,
            x: xVals,
            y: yVals,
            type: 'surface',
            name: `λ=${field.label(c)}`,
            showscale: cFilter !== null || c === 0,
            colorscale: colorscales[c % colorscales.length],
            opacity: 0.85,
            contours: {
                z: { show: true, usecolormap: true, highlightcolor: "limegreen", project: { z: true } }
            }
        });
    }

    const layout = {
        margin: { l: 0, r: 0, b: 0, t: 0 },
        paper_bgcolor: 'transparent',
        plot_bgcolor: 'transparent',
        scene: {
            camera: {
                up: { x: 0, y: 1, z: 0 },
                eye: { x: 1.5, y: 1.2, z: 1.5 }
            },
            aspectmode: 'cube',
            xaxis: {
                title: { text: 'X (x)', font: { color: axisColor, size: 16, family: 'sans-serif' } },
                tickmode: 'array',
                tickvals: xVals,
                ticktext: tickTexts,
                tickfont: { color: textColor, size: 13 },
                linecolor: axisColor,
                tickcolor: axisColor,
                zeroline: false,
                showline: true,
                showgrid: true,
                gridcolor: isLightMode ? 'rgba(0,0,0,0.08)' : 'rgba(255,255,255,0.08)',
                showbackground: false,
                ticks: 'outside',
                ticklen: 5,
                tickwidth: 2,
                showspikes: false
            },
            yaxis: {
                title: { text: 'Y (y)', font: { color: axisColorG, size: 16, family: 'sans-serif' } },
                tickmode: 'array',
                tickvals: yVals,
                ticktext: tickTexts,
                tickfont: { color: textColorG, size: 13 },
                linecolor: axisColorG,
                tickcolor: axisColorG,
                zeroline: false,
                showline: true,
                showgrid: true,
                gridcolor: isLightMode ? 'rgba(0,0,0,0.08)' : 'rgba(255,255,255,0.08)',
                showbackground: false,
                ticks: 'outside',
                ticklen: 5,
                tickwidth: 2,
                showspikes: false
            },
            zaxis: {
                title: { text: 'Z (z)', font: { color: axisColorB, size: 16, family: 'sans-serif' } },
                tickmode: 'array',
                tickvals: xVals,
                ticktext: tickTexts,
                tickfont: { color: textColorB, size: 13 },
                linecolor: axisColorB,
                tickcolor: axisColorB,
                zeroline: false,
                showline: true,
                showgrid: true,
                gridcolor: isLightMode ? 'rgba(0,0,0,0.08)' : 'rgba(255,255,255,0.08)',
                showbackground: false,
                ticks: 'outside',
                ticklen: 5,
                tickwidth: 2,
                showspikes: false
            },
            annotations: [],
            bgcolor: 'transparent'
        }
    };

    Plotly.purge(g);
    Plotly.newPlot(g, data, layout, { responsive: true, displayModeBar: false }).then(() => {
        if (g && window.Plotly && Plotly.Plots && Plotly.Plots.resize) {
            Plotly.Plots.resize(g);
        }
        // Attach listeners immediately for real-time camera sync
        if (g._plotlyRelayoutHandler) {
            try { g.removeListener('plotly_relayout', g._plotlyRelayoutHandler); } catch (e) { }
            try { g.removeListener('plotly_relayouting', g._plotlyRelayoutHandler); } catch (e) { }
        }
        g._plotlyRelayoutHandler = (eventData) => broadcastFromPlotly(g, eventData);
        g.on('plotly_relayout', g._plotlyRelayoutHandler);
        // plotly_relayouting fires continuously during drag (Plotly 2.x)
        try { g.on('plotly_relayouting', g._plotlyRelayoutHandler); } catch (e) { }

        // Fallback: poll _fullLayout.scene.camera during pointer drag
        if (g._plotlyPointerMoveHandler) {
            g.removeEventListener('pointermove', g._plotlyPointerMoveHandler);
        }
        let _dragging = false;
        g.addEventListener('pointerdown', () => { _dragging = true; });
        g.addEventListener('pointerup', () => { _dragging = false; });
        g._plotlyPointerMoveHandler = () => {
            if (!_dragging || !syncCameras) return;
            broadcastFromPlotly(g, {});
        };
        g.addEventListener('pointermove', g._plotlyPointerMoveHandler);
    });
}

function render2D(targetSvgEl, tabCurveStyle) {
    targetSvgEl = targetSvgEl || svgEl;
    tabCurveStyle = tabCurveStyle || curveStyle;
    const field = currentField;
    const q = field.q;

    const coefA = currentA;
    const coefB = currentB;
    eqDisplay2D.textContent = `2D: λ = ${field.label(coefA)}·x + ${field.label(coefB)}·y`;

    const hostWidth = targetSvgEl.parentElement ? targetSvgEl.parentElement.clientWidth : 0;
    const maxPlot = hostWidth > 0
        ? Math.max(380, Math.min(620, hostWidth - 150))
        : 380;
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

    targetSvgEl.setAttribute('viewBox', `0 0 ${svgW} ${svgH}`);
    targetSvgEl.innerHTML = '';

    const px = originMode === 'desc'
        ? xi => offX + (numGaps - xi) * gridSize
        : xi => offX + xi * gridSize;
    const py = originMode === 'desc'
        ? yi => offY + yi * gridSize
        : yi => offY + (numGaps - yi) * gridSize;

    for (let i = 0; i < q; i++) {
        const lbl = field.label(i);
        targetSvgEl.appendChild(makeText(lbl, px(i), offY + numGaps * gridSize + fontSize + 8, 'axis-label', fontSize));
        const ty = makeText(lbl, offX - 26, py(i) + fontSize * 0.38, 'axis-label', fontSize);
        ty.setAttribute('text-anchor', 'end');
        targetSvgEl.appendChild(ty);
    }

    const cFilter = selectedValue === 'all' ? null : parseInt(selectedValue, 10);

    if (coefA !== 0 || coefB !== 0) {
        for (let k = 0; k < q; k++) {
            if (cFilter !== null && k !== cFilter) continue;

            const pts = collectAffineLinePoints2D(field, coefA, coefB, k);
            if (pts.length === 0) continue;

            const color = PALETTE[k % PALETTE.length];
            const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
            g.setAttribute('class', 'line-group');

            const isStraight = coefA === 0 || coefB === 0;
            const d = isStraight
                ? (() => {
                    let s = `M ${px(pts[0].x)} ${py(pts[0].y)}`;
                    for (let i = 1; i < pts.length; i++) s += ` L ${px(pts[i].x)} ${py(pts[i].y)}`;
                    return s;
                })()
                : buildSVGPath(pts, px, py, tabCurveStyle);

            const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
            path.setAttribute('d', d);
            path.setAttribute('class', 'line-path');
            path.setAttribute('stroke', color);
            path.setAttribute('stroke-width', strokeW);
            path.classList.add('redraw');
            g.appendChild(path);

            const last = pts[pts.length - 1];
            const lblStr = `λ=${field.label(k)}`;
            const lbl = makeText(lblStr, px(last.x) + 7, py(last.y) + fontSize * 0.3, 'equation-label', Math.max(8, fontSize - 1));
            lbl.setAttribute('fill', color);
            lbl.setAttribute('text-anchor', 'start');
            g.appendChild(lbl);

            targetSvgEl.appendChild(g);
        }
    }

    if (cFilter === null) {
        for (let y = 0; y < q; y++) {
            for (let x = 0; x < q; x++) {
                const dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                dot.setAttribute('cx', px(x));
                dot.setAttribute('cy', py(y));
                dot.setAttribute('r', dotR);
                dot.setAttribute('class', 'grid-point');
                targetSvgEl.appendChild(dot);
            }
        }
    } else {
        for (let y = 0; y < q; y++) {
            for (let x = 0; x < q; x++) {
                const lam = field.add(field.mul(coefA, x), field.mul(coefB, y));
                if (lam !== cFilter) continue;
                const dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                dot.setAttribute('cx', px(x));
                dot.setAttribute('cy', py(y));
                dot.setAttribute('r', dotR);
                dot.setAttribute('class', 'grid-point');
                dot.setAttribute('fill', PALETTE[cFilter % PALETTE.length]);
                targetSvgEl.appendChild(dot);
            }
        }
    }
}

function renderLatinSquare(field, q) {
    latinSquareEl.innerHTML = '';

    const maxMatrixSize = Math.min(520, Math.max(320, q * 56));
    const gap = Math.max(3, Math.min(8, Math.floor(maxMatrixSize / q / 10)));
    const cellSize = Math.max(28, Math.min(90, Math.floor((maxMatrixSize - (q - 1) * gap) / q)));
    const cellFont = Math.max(8, Math.min(16, Math.floor(cellSize * 0.4)));

    latinSquareEl.style.gap = `${gap}px`;
    latinSquareEl.style.gridTemplateColumns = `repeat(${q}, ${cellSize}px)`;
    latinSquareEl.style.width = `${q * cellSize + (q - 1) * gap + 28}px`;
    latinSquareEl.style.height = `${q * cellSize + (q - 1) * gap + 28}px`;

    const cFilter = selectedValue === 'all' ? null : parseInt(selectedValue, 10);
    const coefA = currentA;
    const coefB = currentB;

    // originMode affects both x-axis and y-axis order:
    const xIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => q - 1 - i)
        : Array.from({ length: q }, (_, i) => i);
    const yIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => i)
        : Array.from({ length: q }, (_, i) => q - 1 - i);

    for (let row = 0; row < q; row++) {
        const yi = yIndices[row];
        for (let col = 0; col < q; col++) {
            const x = xIndices[col];
            const lam = field.add(field.mul(coefA, x), field.mul(coefB, yi));

            const cellText = (cFilter !== null && lam !== cFilter) ? '' : field.label(lam);
            const detailText = cellText
                ? `x = ${field.label(x)} | y = ${field.label(yi)} | λ = ${field.label(lam)}`
                : '';
            latinSquareEl.appendChild(makeCell(cellText, cellText ? PALETTE[lam % PALETTE.length] : '#333', cellSize, cellFont, detailText));
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

function makeCell(text, color, size, fontSize, tooltipText = '') {
    const div = document.createElement('div');
    div.className = 'latin-cell';
    div.textContent = text;
    div.style.color = color;
    div.style.borderColor = color + '33';
    div.style.width = size + 'px';
    div.style.height = size + 'px';
    div.style.fontSize = fontSize + 'px';

    if (tooltipText && latinTooltip) {
        div.addEventListener('mouseenter', () => {
            latinTooltip.textContent = tooltipText;
            latinTooltip.classList.add('visible');
        });
        div.addEventListener('mousemove', (event) => {
            const offsetX = 0;
            const offsetY = 24;
            latinTooltip.style.left = `${event.clientX + offsetX}px`;
            latinTooltip.style.top = `${event.clientY - offsetY}px`;
        });
        div.addEventListener('mouseleave', () => {
            latinTooltip.classList.remove('visible');
        });
    }

    return div;
}

// ===================================================================
// Browser-tab 3D panel system
// ===================================================================

function initBrowserTabs() {
    const addMenu = document.getElementById('browser-tab-add-menu');

    addMenu.querySelectorAll('.browser-tab-add-menu-item').forEach(item => {
        item.addEventListener('click', (e) => {
            e.stopPropagation();
            const type = item.dataset.panel;
            if (openTabs.length < MAX_PANELS) {
                addBrowserTab(type);
            }
        });
    });
}

function refreshAddMenu() {
    const addMenu = document.getElementById('browser-tab-add-menu');
    const count = openTabs.length;
    addMenu.querySelectorAll('.browser-tab-add-menu-item').forEach(item => {
        item.classList.toggle('disabled', count >= MAX_PANELS);
    });
}

function addBrowserTab(type) {
    if (openTabs.length >= MAX_PANELS) return;
    const id = createTabId();

    // Create panel DOM
    const def = PANEL_DEFS[type];
    const panelEl = document.createElement('div');
    panelEl.className = 'browser-panel-item browser-3d-panel-item';
    panelEl.id = `browser-panel-${id}`;
    if (type === 'linesPlanes') panelEl.classList.add('browser-panel-lines-planes');

    const titleEl = document.createElement('div');
    titleEl.className = 'browser-panel-title';
    titleEl.textContent = def.title;
    panelEl.appendChild(titleEl);

    let visualizer = null;
    let plotlyEl = null;
    let showNum = false;
    let tabCurveStyle = curveStyle;  // per-tab independent curve style
    let latinPanelEl = null;
    let sliceLabelEl = null;

    if (type === 'plotly') {
        plotlyEl = document.createElement('div');
        plotlyEl.style.cssText = 'width:100%;height:100%;border-radius:10px;overflow:hidden;background:transparent;';
        panelEl.appendChild(plotlyEl);
    } else {
        const innerEl = document.createElement('div');
        innerEl.className = 'browser-panel-inner';
        innerEl.style.cssText = 'width:100%;height:100%;min-height:260px;';

        panelEl.appendChild(innerEl);

        // Show-number toggle for latin only
        if (type === 'latin') {
            const toggleEl = document.createElement('div');
            toggleEl.className = 'browser-panel-show-num';
            toggleEl.innerHTML = `<input type="checkbox" id="show-num-${id}"><label for="show-num-${id}">顯示數字</label>`;
            panelEl.appendChild(toggleEl);
            toggleEl.querySelector('input').addEventListener('change', (e) => {
                const tab = openTabs.find(t => t.id === id);
                if (tab) {
                    tab.showNum = e.target.checked;
                    renderSingleTab(tab);
                }
            });
        }

        // Curve style overlay for lines / 線圖＋疊層
        if (type === 'lines' || type === 'linesPlanes') {
            const curveOverlay = document.createElement('div');
            curveOverlay.className = 'panel-curve-style-overlay';
            curveOverlay.innerHTML = `
                <div class="panel-curve-style-title" title="線條畫法">≡</div>
                <button class="panel-curve-btn${curveStyle === 'bezier' ? ' active' : ''}" data-style="bezier">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 C10,18 14,2 20,2 C26,2 30,18 38,18" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>Bézier
                </button>
                <button class="panel-curve-btn${curveStyle === 'smoothstep' ? ' active' : ''}" data-style="smoothstep">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 C12,18 16,2 38,2" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>Smoothstep
                </button>
                <button class="panel-curve-btn${curveStyle === 'catmull' ? ' active' : ''}" data-style="catmull">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,14 C8,18 12,2 20,10 C28,18 32,2 38,6" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>Catmull-Rom
                </button>
                <button class="panel-curve-btn${curveStyle === 'arc' ? ' active' : ''}" data-style="arc">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 A18 18 0 0 1 20,2 A18 18 0 0 1 38,18" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>弧線
                </button>
                <button class="panel-curve-btn${curveStyle === 'step' ? ' active' : ''}" data-style="step">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><polyline points="2,18 14,18 14,8 26,8 26,2 38,2" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>階梯
                </button>
                <button class="panel-curve-btn${curveStyle === 'straight' ? ' active' : ''}" data-style="straight">
                    <span class="curve-icon"><svg viewBox="0 0 40 20" width="34" height="17"><polyline points="2,18 14,2 26,14 38,6" stroke="currentColor" stroke-width="2" fill="none"/></svg></span>直線
                </button>
            `;
            curveOverlay.querySelectorAll('.panel-curve-btn').forEach(btn => {
                btn.addEventListener('click', () => {
                    const style = btn.dataset.style;
                    // Update only this tab's curveStyle
                    const tab = openTabs.find(t => t.id === id);
                    if (tab) {
                        tab.curveStyle = style;
                        // Update button states within this overlay only
                        curveOverlay.querySelectorAll('.panel-curve-btn').forEach(b => b.classList.toggle('active', b.dataset.style === style));
                        renderSingleTab(tab);
                    }
                });
            });
            panelEl.appendChild(curveOverlay);
        }

        // Create Three.js visualizer after appending to DOM
        document.getElementById('browser-panels-row').appendChild(panelEl);
        visualizer = new Visualizer3D(innerEl);
    }

    if (type !== 'plotly') {
        // already appended above for Three.js init
    } else {
        document.getElementById('browser-panels-row').appendChild(panelEl);
    }

    const tabEntry = {
        id,
        type,
        visualizer,
        plotlyEl,
        showNum,
        panelEl,
        curveStyle: tabCurveStyle,
        latinPanelEl,
        sliceLabelEl,
    };
    openTabs.push(tabEntry);

    // 同步新 panel 相機到現有第一個 Three.js panel 的視角
    if (visualizer) {
        const srcTab = openTabs.slice(0, -1).find(t => t.visualizer && t.type !== 'plotly');
        if (srcTab) {
            visualizer.camera.position.copy(srcTab.visualizer.camera.position);
            visualizer.camera.quaternion.copy(srcTab.visualizer.camera.quaternion);
            visualizer.controls.target.copy(srcTab.visualizer.controls.target);
            visualizer.controls.update();
            // 預先標記 lastQ，避免 render() 覆寫已同步的相機位置
            if (currentField) {
                visualizer.lastQ = currentField.q;
                visualizer._lastOriginMode = originMode;
            }
        }
    }

    renderTabBar();
    update3DPanelLayout();
    if (currentField) renderSingleTab(tabEntry);

    // Trigger resize
    setTimeout(() => window.dispatchEvent(new Event('resize')), 50);
}

function removeBrowserTab(id) {
    const idx = openTabs.findIndex(t => t.id === id);
    if (idx === -1) return;
    const tab = openTabs[idx];

    // Dispose Three.js if applicable
    if (tab.visualizer) {
        if (tab.visualizer._ro) tab.visualizer._ro.disconnect();
        tab.visualizer.clear();
        if (tab.visualizer.renderer) {
            tab.visualizer.renderer.dispose();
            tab.visualizer.renderer.domElement.remove();
        }
    }
    if (tab.plotlyEl && window.Plotly) {
        Plotly.purge(tab.plotlyEl);
    }
    tab.panelEl.remove();
    openTabs.splice(idx, 1);
    renderTabBar();
    update3DPanelLayout();
}

function renderTabBar() {
    const list = document.getElementById('browser-tabs-list');
    list.innerHTML = '';

    openTabs.forEach(tab => {
        const def = PANEL_DEFS[tab.type];
        const tabEl = document.createElement('div');
        tabEl.className = 'browser-tab active'; // all visible tabs are "active" (always showing)
        tabEl.dataset.id = tab.id;

        const labelEl = document.createElement('span');
        labelEl.className = 'browser-tab-label';
        labelEl.textContent = def.label;

        const closeEl = document.createElement('span');
        closeEl.className = 'browser-tab-close';
        closeEl.textContent = '✕';
        closeEl.title = '關閉';
        closeEl.addEventListener('click', (e) => {
            e.stopPropagation();
            removeBrowserTab(tab.id);
        });

        tabEl.appendChild(labelEl);
        tabEl.appendChild(closeEl);
        list.appendChild(tabEl);
    });

    // Update add button disabled state
    const addBtn = document.getElementById('browser-tab-add');
    const addWrapper = addBtn && addBtn.closest('.browser-tab-add-wrapper');
    const full = openTabs.length >= MAX_PANELS;
    addBtn.style.opacity = full ? '0.3' : '1';
    addBtn.style.pointerEvents = full ? 'none' : 'auto';
    if (addWrapper) addWrapper.style.pointerEvents = full ? 'none' : 'auto';
    // Also update disabled class on menu items
    const addMenu = document.getElementById('browser-tab-add-menu');
    if (addMenu) {
        addMenu.querySelectorAll('.browser-tab-add-menu-item').forEach(item => {
            item.classList.toggle('disabled', full);
        });
    }
}

function renderSingleTab(tab) {
    if (!currentField) return;
    const axisFilters = {
        x: selectedAxisX,
        y: selectedAxisY,
        z: selectedAxisZ,
        c: selectedValue
    };
    if (tab.type === 'plotly') {
        renderPlotlySurface(currentField, currentA, currentB, currentC, axisFilters, tab.plotlyEl);
    } else {
        const mode = PANEL_DEFS[tab.type].mode;
        tab.visualizer.render(mode, currentField, tab.showNum, currentA, currentB, currentC, axisFilters, tab.curveStyle || curveStyle);



        // Fade-in
        const canvas = tab.visualizer.renderer && tab.visualizer.renderer.domElement;
        if (canvas) {
            canvas.style.opacity = '0';
            setTimeout(() => { canvas.style.opacity = '1'; }, 40);
        }
    }
}

function renderAllOpenTabs() {
    openTabs.forEach(tab => renderSingleTab(tab));
}

function update3DPanelLayout() {
    const rowEl = document.getElementById('browser-panels-row');
    if (!rowEl) return;
    const count = openTabs.length;
    if (count <= 0) return;

    const gap = 12;
    rowEl.style.justifyContent = count === 1 ? 'center' : 'space-between';

    openTabs.forEach(tab => {
        if (!tab.panelEl) return;
        if (count === 1) {
            tab.panelEl.style.flex = '0 1 min(620px, 82%)';
        } else {
            tab.panelEl.style.flex = `1 1 calc((100% - ${(count - 1) * gap}px) / ${count})`;
        }
    });
}

// ===== 相機同步 =====

// Three.js camera position → Plotly eye (normalized to distance ~2)
function threeToPlotlyEye(cam, target) {
    const dx = cam.position.x - target.x;
    const dy = cam.position.y - target.y;
    const dz = cam.position.z - target.z;
    const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    const scale = dist > 0 ? 2.0 / dist : 1;
    return { x: dx * scale, y: dy * scale, z: dz * scale };
}

// Plotly eye → Three.js camera offset from target
function plotlyEyeToThree(eye, target, desiredDist) {
    const len = Math.sqrt(eye.x * eye.x + eye.y * eye.y + eye.z * eye.z) || 1;
    const scale = desiredDist / len;
    return new THREE.Vector3(
        target.x + eye.x * scale,
        target.y + eye.y * scale,
        target.z + eye.z * scale
    );
}

let _syncLock = false; // prevent feedback loops

function broadcastCamera(source) {
    if (_syncLock) return;
    _syncLock = true;

    const eye = threeToPlotlyEye(source.camera, source.controls.target);
    const dist = source.camera.position.distanceTo(source.controls.target);

    openTabs.forEach(tab => {
        if (tab.type === 'plotly') {
            // Sync to Plotly via relayout
            if (tab.plotlyEl && tab.plotlyEl._fullLayout) {
                Plotly.relayout(tab.plotlyEl, { 'scene.camera': { eye, up: { x: 0, y: 1, z: 0 }, center: { x: 0, y: 0, z: 0 } } });
            }
        } else {
            if (!tab.visualizer || tab.visualizer === source) return;
            const dst = tab.visualizer;
            dst.camera.position.copy(source.camera.position);
            dst.camera.quaternion.copy(source.camera.quaternion);
            dst.controls.target.copy(source.controls.target);
            dst.controls.update();
        }
    });

    _syncLock = false;
}

function broadcastFromPlotly(plotlyEl, eventData) {
    if (_syncLock) return;
    if (!syncCameras) return;

    // During drag, _fullLayout always has the most up-to-date camera
    let eye = null;
    const layout = plotlyEl._fullLayout;
    if (layout && layout.scene && layout.scene.camera && layout.scene.camera.eye) {
        eye = layout.scene.camera.eye;
    } else if (eventData && eventData['scene.camera'] && eventData['scene.camera'].eye) {
        eye = eventData['scene.camera'].eye;
    } else if (eventData && eventData['scene.camera.eye']) {
        eye = eventData['scene.camera.eye'];
    }
    if (!eye || (eye.x === undefined)) return;

    _syncLock = true;

    openTabs.forEach(tab => {
        if (tab.type === 'plotly') {
            if (tab.plotlyEl === plotlyEl) return; // skip self
            if (tab.plotlyEl && tab.plotlyEl._fullLayout) {
                Plotly.relayout(tab.plotlyEl, {
                    'scene.camera': { eye, up: { x: 0, y: 1, z: 0 }, center: { x: 0, y: 0, z: 0 } }
                });
            }
        } else {
            if (!tab.visualizer) return;
            const dst = tab.visualizer;
            const desiredDist = dst.camera.position.distanceTo(dst.controls.target) || 20;
            const center = dst.controls.target;
            dst.camera.position.copy(plotlyEyeToThree(eye, center, desiredDist));
            dst.camera.lookAt(center);
            dst.controls.update();
        }
    });

    _syncLock = false;
}

// ===================================================================
// Browser-tab 2D panel system
// ===================================================================

function initBrowserTabs2D() {
    const addMenu = document.getElementById('browser-2d-tab-add-menu');
    if (!addMenu) return;
    addMenu.querySelectorAll('.browser-2d-tab-add-menu-item').forEach(item => {
        item.addEventListener('click', (e) => {
            e.stopPropagation();
            const type = item.dataset.panel;
            if (open2DTabs.length < MAX_2D_PANELS) {
                add2DTab(type);
            }
        });
    });
    // close menu on outside click
    document.addEventListener('click', (e) => {
        const wrapper = document.getElementById('browser-2d-tab-add-wrapper');
        if (wrapper && !wrapper.contains(e.target)) {
            addMenu.classList.remove('open');
        }
    });
    const addBtn = document.getElementById('browser-2d-tab-add');
    if (addBtn) {
        addBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            addMenu.classList.toggle('open');
        });
    }
}

function add2DTab(type) {
    if (open2DTabs.length >= MAX_2D_PANELS) return;
    const id = create2DTabId();
    const def = PANEL_2D_DEFS[type];

    const panelEl = document.createElement('div');
    panelEl.className = 'browser-panel-item browser-2d-panel-item';
    panelEl.id = `browser-2d-panel-${id}`;

    const titleEl = document.createElement('div');
    titleEl.className = 'browser-panel-title';
    titleEl.textContent = def.title;
    panelEl.appendChild(titleEl);

    let tabCurveStyle = curveStyle;

    if (type === 'line') {
        // SVG line chart with its own curve-style overlay
        const svgWrapper = document.createElement('div');
        svgWrapper.className = 'viz-2d-wrapper';
        svgWrapper.style.cssText = 'width:100%;flex:1;display:flex;align-items:center;justify-content:center;position:relative;';

        const newSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        newSvg.setAttribute('overflow', 'visible');
        newSvg.style.cssText = 'flex:1;width:100%;max-width:100%;';
        svgWrapper.appendChild(newSvg);

        // Curve style overlay
        const overlay = document.createElement('div');
        overlay.className = 'panel-curve-style-overlay';
        overlay.innerHTML = buildCurveStyleOverlayHTML(tabCurveStyle);
        overlay.querySelectorAll('.panel-curve-btn').forEach(btn => {
            btn.addEventListener('click', () => {
                const style = btn.dataset.style;
                const tab = open2DTabs.find(t => t.id === id);
                if (tab) {
                    tab.curveStyle = style;
                    overlay.querySelectorAll('.panel-curve-btn').forEach(b => b.classList.toggle('active', b.dataset.style === style));
                    render2D(tab.svgEl, tab.curveStyle);
                }
            });
        });
        svgWrapper.appendChild(overlay);
        panelEl.appendChild(svgWrapper);

        document.getElementById('browser-2d-panels-row').appendChild(panelEl);
        const tabEntry = { id, type, svgEl: newSvg, latinEl: null, panelEl, curveStyle: tabCurveStyle };
        open2DTabs.push(tabEntry);
        render2DTab(tabEntry);

    } else if (type === 'latin') {
        const latinContainer = document.createElement('div');
        latinContainer.className = 'visualization latin-square-container';
        latinContainer.style.cssText = 'flex:1;display:flex;flex-direction:column;align-items:center;justify-content:center;padding:16px;';

        const latinTitle = document.createElement('div');
        latinTitle.className = 'latin-square-title';
        latinTitle.style.cssText = 'margin-bottom:12px;font-size:1.05em;color:var(--primary-color);';
        latinTitle.textContent = '拉丁方陣 (2D: c = y - ax)';
        latinContainer.appendChild(latinTitle);

        const latinGrid = document.createElement('div');
        latinGrid.className = 'latin-square';
        latinContainer.appendChild(latinGrid);
        panelEl.appendChild(latinContainer);

        document.getElementById('browser-2d-panels-row').appendChild(panelEl);
        const tabEntry = { id, type, svgEl: null, latinEl: latinGrid, panelEl, curveStyle: tabCurveStyle };
        open2DTabs.push(tabEntry);
        render2DTab(tabEntry);
    }

    render2DTabBar();
    update2DPanelLayout();
    setTimeout(() => window.dispatchEvent(new Event('resize')), 50);
}

function buildCurveStyleOverlayHTML(activeStyle) {
    const styles = [
        { key: 'bezier', label: 'Bézier', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 C10,18 14,2 20,2 C26,2 30,18 38,18" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
        { key: 'smoothstep', label: 'Smoothstep', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 C12,18 16,2 38,2" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
        { key: 'catmull', label: 'Catmull-Rom', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,14 C8,18 12,2 20,10 C28,18 32,2 38,6" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
        { key: 'arc', label: '弧線', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><path d="M2,18 A18 18 0 0 1 20,2 A18 18 0 0 1 38,18" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
        { key: 'step', label: '階梯', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><polyline points="2,18 14,18 14,8 26,8 26,2 38,2" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
        { key: 'straight', label: '直線', icon: `<svg viewBox="0 0 40 20" width="34" height="17"><polyline points="2,18 14,2 26,14 38,6" stroke="currentColor" stroke-width="2" fill="none"/></svg>` },
    ];
    return `<div class="panel-curve-style-title" title="線條畫法">≡</div>` +
        styles.map(s => `<button class="panel-curve-btn${s.key === activeStyle ? ' active' : ''}" data-style="${s.key}">
            <span class="curve-icon">${s.icon}</span>${s.label}</button>`).join('');
}

function remove2DTab(id) {
    const idx = open2DTabs.findIndex(t => t.id === id);
    if (idx === -1) return;
    const tab = open2DTabs[idx];
    tab.panelEl.remove();
    open2DTabs.splice(idx, 1);
    render2DTabBar();
    update2DPanelLayout();
}

function render2DTabBar() {
    const list = document.getElementById('browser-2d-tabs-list');
    if (!list) return;
    list.innerHTML = '';
    open2DTabs.forEach(tab => {
        const def = PANEL_2D_DEFS[tab.type];
        const tabEl = document.createElement('div');
        tabEl.className = 'browser-tab active';
        tabEl.dataset.id = tab.id;

        const labelEl = document.createElement('span');
        labelEl.className = 'browser-tab-label';
        labelEl.textContent = def.label;

        const closeEl = document.createElement('span');
        closeEl.className = 'browser-tab-close';
        closeEl.textContent = '✕';
        closeEl.title = '關閉';
        closeEl.addEventListener('click', (e) => {
            e.stopPropagation();
            remove2DTab(tab.id);
        });
        tabEl.appendChild(labelEl);
        tabEl.appendChild(closeEl);
        list.appendChild(tabEl);
    });

    const addBtn = document.getElementById('browser-2d-tab-add');
    const full = open2DTabs.length >= MAX_2D_PANELS;
    if (addBtn) {
        addBtn.style.opacity = full ? '0.3' : '1';
        addBtn.style.pointerEvents = full ? 'none' : 'auto';
    }
    const addMenu = document.getElementById('browser-2d-tab-add-menu');
    if (addMenu) {
        addMenu.querySelectorAll('.browser-2d-tab-add-menu-item').forEach(item => {
            item.classList.toggle('disabled', full);
        });
    }
}

function update2DPanelLayout() {
    const rowEl = document.getElementById('browser-2d-panels-row');
    if (!rowEl) return;
    const count = open2DTabs.length;
    if (count <= 0) return;

    const gap = 12;
    rowEl.style.justifyContent = 'center';

    open2DTabs.forEach(tab => {
        if (!tab.panelEl) return;
        if (count === 1) {
            tab.panelEl.style.flex = '0 1 min(760px, 90%)';
            tab.panelEl.style.maxWidth = '760px';
        } else {
            tab.panelEl.style.flex = `0 1 calc((100% - ${(count - 1) * gap}px) / ${count})`;
            tab.panelEl.style.maxWidth = count === 2 ? '620px' : 'none';
        }
    });
}

function render2DTab(tab) {
    if (!currentField) return;
    if (tab.type === 'line' && tab.svgEl) {
        render2D(tab.svgEl, tab.curveStyle);
    } else if (tab.type === 'latin' && tab.latinEl) {
        renderLatinSquareToEl(tab.latinEl, currentField, currentField.q);
    }
}

function renderAllOpen2DTabs() {
    open2DTabs.forEach(tab => render2DTab(tab));
}

function renderLatinSquareToEl(targetEl, field, q) {
    targetEl.innerHTML = '';
    const maxMatrixSize = Math.min(480, Math.max(300, q * 52));
    const gap = Math.max(3, Math.min(8, Math.floor(maxMatrixSize / q / 10)));
    const cellSize = Math.max(26, Math.min(80, Math.floor((maxMatrixSize - (q - 1) * gap) / q)));
    const cellFont = Math.max(8, Math.min(15, Math.floor(cellSize * 0.4)));

    targetEl.style.gap = `${gap}px`;
    targetEl.style.gridTemplateColumns = `repeat(${q}, ${cellSize}px)`;
    targetEl.style.width = `${q * cellSize + (q - 1) * gap + 28}px`;
    targetEl.style.height = `${q * cellSize + (q - 1) * gap + 28}px`;

    const cFilter = selectedValue === 'all' ? null : parseInt(selectedValue, 10);
    const coefA = currentA;
    const coefB = currentB;

    const xIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => q - 1 - i)
        : Array.from({ length: q }, (_, i) => i);
    const yIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => i)
        : Array.from({ length: q }, (_, i) => q - 1 - i);

    for (let row = 0; row < q; row++) {
        const yi = yIndices[row];
        for (let col = 0; col < q; col++) {
            const x = xIndices[col];
            const lam = field.add(field.mul(coefA, x), field.mul(coefB, yi));
            const cellText = (cFilter !== null && lam !== cFilter) ? '' : field.label(lam);
            const detailText = cellText ? `x = ${field.label(x)} | y = ${field.label(yi)} | λ = ${field.label(lam)}` : '';
            targetEl.appendChild(makeCell(cellText, cellText ? PALETTE[lam % PALETTE.length] : '#333', cellSize, cellFont, detailText));
        }
    }
}

/** 固定 X 時，(Y,Z) 網格上的 λ = ax+by+cz（橫向 Y、縱向 Z，排列習慣同 2D Latin） */
function render3DSliceMatrixPanel(targetEl, field, q, xSlice, coefA, coefB, coefC, cFilter) {
    targetEl.innerHTML = '';
    const maxMatrixSize = Math.min(440, Math.max(280, q * 48));
    const gap = Math.max(3, Math.min(8, Math.floor(maxMatrixSize / q / 10)));
    const cellSize = Math.max(24, Math.min(76, Math.floor((maxMatrixSize - (q - 1) * gap) / q)));
    const cellFont = Math.max(8, Math.min(14, Math.floor(cellSize * 0.38)));

    targetEl.style.gap = `${gap}px`;
    targetEl.style.gridTemplateColumns = `repeat(${q}, ${cellSize}px)`;
    targetEl.style.width = `${q * cellSize + (q - 1) * gap + 24}px`;
    targetEl.style.height = `${q * cellSize + (q - 1) * gap + 24}px`;

    const yIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => q - 1 - i)
        : Array.from({ length: q }, (_, i) => i);
    const zIndices = originMode === 'desc'
        ? Array.from({ length: q }, (_, i) => i)
        : Array.from({ length: q }, (_, i) => q - 1 - i);

    for (let row = 0; row < q; row++) {
        const zi = zIndices[row];
        for (let col = 0; col < q; col++) {
            const yj = yIndices[col];
            const lam = affineLambda(field, coefA, coefB, coefC, xSlice, yj, zi);
            const cellText = (cFilter !== null && lam !== cFilter) ? '' : field.label(lam);
            const detailText = cellText
                ? `x = ${field.label(xSlice)} | y = ${field.label(yj)} | z = ${field.label(zi)} | λ = ${field.label(lam)}`
                : '';
            targetEl.appendChild(makeCell(cellText, cellText ? PALETTE[lam % PALETTE.length] : '#333', cellSize, cellFont, detailText));
        }
    }
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
