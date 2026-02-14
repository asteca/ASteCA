// This script is used to load the Umami analytics script asynchronously.
const script = document.createElement('script');
script.defer = true;
script.src = 'https://cloud.umami.is/script.js';
script.setAttribute('data-website-id', '8cd88b63-caf5-4d82-98e2-8f765e390008');
document.head.appendChild(script);
