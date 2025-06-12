<script setup>
import { ref, watch, onMounted } from "vue";

defineProps({
    isSidebarOpen: {
        type: Boolean,
        required: true,
    },
});

defineEmits(["toggle-sidebar", "goto-home"]);

const sidebarWidth = ref(25);
const isResizing = ref(false);
const sidebarRef = ref(null);
const isDarkMode = ref(
    localStorage.getItem("darkMode") === "true" ||
        (window.matchMedia &&
            window.matchMedia("(prefers-color-scheme: dark)").matches),
);

// Function to toggle theme
function toggleTheme() {
    isDarkMode.value = !isDarkMode.value;
    localStorage.setItem("theme", isDarkMode.value ? "dark" : "light");
    document.documentElement.setAttribute(
        "data-theme",
        isDarkMode.value ? "dark" : "light",
    );
}

const handleMouseDown = (event) => {
    isResizing.value = true;
    document.addEventListener("mousemove", handleMouseMove);
    document.addEventListener("mouseup", handleMouseUp);
};

const handleMouseMove = (event) => {
    if (!isResizing.value) return;
    const newWidth = (event.clientX / window.innerWidth) * 100;
    sidebarWidth.value = Math.min(Math.max(newWidth, 20), 35);
};

const handleMouseUp = () => {
    isResizing.value = false;
    document.removeEventListener("mousemove", handleMouseMove);
    document.removeEventListener("mouseup", handleMouseUp);
};

watch(
    () => isResizing.value,
    (newVal) => {
        if (sidebarRef.value) {
            if (newVal) {
                sidebarRef.value.classList.add("resizing");
            } else {
                sidebarRef.value.classList.remove("resizing");
            }
        }
    },
);

onMounted(() => {
    const initialTheme = localStorage.getItem("theme");
    isDarkMode.value = initialTheme === "dark";
    document.documentElement.setAttribute(
        "data-theme",
        isDarkMode.value ? "dark" : "light",
    );
    if (sidebarRef.value) {
        sidebarRef.value.style.width = `${sidebarWidth.value}%`;
    }
});

watch(
    () => sidebarWidth.value,
    (newVal) => {
        if (sidebarRef.value) {
            sidebarRef.value.style.width = `${newVal}%`;
        }
    },
);

// Watch for changes to isDarkMode to update document class list
watch(
    () => isDarkMode.value,
    (newVal) => {
        if (newVal) {
            document.documentElement.classList.add("dark");
        } else {
            document.documentElement.classList.remove("dark");
        }
    },
    { immediate: true },
);
</script>

<template>
    <aside
        ref="sidebarRef"
        class="my-sidebar"
        :class="{ 'my-sidebar-open': isSidebarOpen }"
    >
        <div class="resize-handle" @mousedown="handleMouseDown"></div>
        <div class="sidebar-header">
            <button class="home-button" @click="$emit('goto-home')">
                <i class="fas fa-home"></i>
            </button>
            <span class="spacer"></span>
            <button
                class="home-button"
                @click="toggleTheme"
                :title="
                    isDarkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'
                "
            >
                <i :class="isDarkMode ? 'fas fa-sun' : 'fas fa-moon'"></i>
            </button>
            <span class="header-spacer"></span>
            <button class="close-button" @click="$emit('toggle-sidebar')">
                <i
                    :class="
                        isSidebarOpen ? 'fas fa-chevron-left' : 'fas fa-bars'
                    "
                ></i>
            </button>
        </div>
        <div class="sidebar-content">
            <slot></slot>
        </div>
    </aside>
</template>

<style scoped>
.resize-handle {
    position: absolute;
    top: 0;
    right: -5px;
    width: 10px;
    height: 100%;
    cursor: col-resize;
    background: transparent;
    z-index: 1001;
}

.my-sidebar {
    position: fixed;
    top: 0;
    left: 0;
    height: 100vh;
    min-width: 20%;
    max-width: 35%;
    background-color: var(--sidebar-background);
    transition:
        transform 0.3s ease-in-out,
        background-color 0.3s ease-in-out,
        box-shadow 0.3s ease-in-out;
    box-shadow: 2px 0 5px var(--shadow-color);
    transform: translateX(-100%);
    display: flex;
    flex-direction: column;
    z-index: 1000;
}

.my-sidebar.resizing {
    transition: none;
}

.my-sidebar-open {
    transform: translateX(0);
}

@media (max-width: 768px) {
    .my-sidebar {
        width: 100% !important;
        max-width: 100%;
        height: 100vh;
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
    }

    .resize-handle {
        display: none;
    }

    .sidebar-header {
        padding: 1rem;
        display: flex;
        align-items: center;
        justify-content: space-between;
        border-bottom: 1px solid #ddd;
    }

    .sidebar-content {
        flex-grow: 1;
        overflow-y: auto;
        padding: 0.5rem 1rem;
    }
}

.sidebar-header {
    padding: 0.5rem;
    display: flex;
    align-items: center;
    justify-content: space-between;
    border-bottom: 1px solid #ddd;
}

.sidebar-content {
    flex-grow: 1;
    overflow-y: auto;
    padding: 0.5rem 1rem;
}

.close-button,
.home-button,
.theme-button {
    background: none;
    border: none;
    padding: 0.25rem;
    cursor: pointer;
    font-size: 1.2em;
    color: var(--text-color);
}

.spacer {
    flex: 1;
}

.header-spacer {
    width: 1.5rem;
}
</style>
