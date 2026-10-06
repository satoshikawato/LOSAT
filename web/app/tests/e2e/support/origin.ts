// The origin of the application that the E2E tests open (`vite preview`): port 4173, or
// LOSAT_WEB_E2E_PORT.
export const E2E_PORT = Number(process.env['LOSAT_WEB_E2E_PORT'] || 4173);
export const E2E_ORIGIN = `http://localhost:${E2E_PORT}`;
