// Emscripten이 생성한 mersenne_wasm.js 불러오기
importScripts('mersenne_wasm.js');

let checkMersenne = null;

// Wasm 모듈 로딩 완료 시점
Module.onRuntimeInitialized = () => {
    // C++의 check_mersenne_prime 함수 바인딩
    checkMersenne = Module.cwrap('check_mersenne_prime', 'number', ['number']);
    self.postMessage({ type: 'READY' });
};

// 메인 스레드(HTML)로부터 연산 요청을 받았을 때
self.onmessage = (e) => {
    const { action, startP, endP } = e.data;

    if (action === 'START' && checkMersenne) {
        for (let p = startP; p <= endP; p++) {
            const isPrime = checkMersenne(p);
            
            if (isPrime === 1) {
                // 메르센 소수 발견!
                self.postMessage({ type: 'FOUND', p: p });
            } else {
                // 탐색 진행 중
                self.postMessage({ type: 'PROGRESS', p: p });
            }
        }
        self.postMessage({ type: 'DONE' });
    }
};
