---
source: /Users/byron/projects/rockbeatspaper/src/routes/index.tsx
relative: rockbeatspaper/src/routes/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute, Link } from "@tanstack/react-router";
import {
  Zap,
  Brain,
  TrendingDown,
  Shield,
  BarChart3,
  Thermometer,
} from "lucide-react";
import { Navbar } from "~/components/Navbar";
import { Button } from "~/components/Button";
import { Card } from "~/components/Card";

export const Route = createFileRoute("/")({
  component: Home,
});

function Home() {
  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-purple-50 to-pink-50 dark:from-gray-900 dark:via-gray-800 dark:to-gray-900">
      <Navbar />

      {/* Hero Section */}
      <section className="relative overflow-hidden">
        <div className="absolute inset-0 bg-gradient-to-r from-blue-600/10 to-purple-600/10 backdrop-blur-3xl"></div>
        <div className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-24 sm:py-32">
          <div className="text-center">
            <h1 className="text-5xl sm:text-6xl lg:text-7xl font-extrabold text-gray-900 dark:text-white mb-6">
              <span className="bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent">
                AI-Powered
              </span>
              <br />
              Thermal Management
            </h1>
            <p className="text-xl sm:text-2xl text-gray-600 dark:text-gray-300 mb-8 max-w-3xl mx-auto">
              Optimize your energy consumption with intelligent thermal analysis.
              Save money, reduce carbon footprint, and maintain perfect comfort.
            </p>
            <div className="flex flex-col sm:flex-row gap-4 justify-center">
              <Link to="/signup">
                <Button size="lg" className="w-full sm:w-auto">
                  Get Started Free
                </Button>
              </Link>
              <Link to="/login">
                <Button variant="outline" size="lg" className="w-full sm:w-auto">
                  Sign In
                </Button>
              </Link>
            </div>
          </div>

          {/* Hero Image */}
          <div className="mt-16 relative">
            <div className="absolute inset-0 bg-gradient-to-r from-blue-600 to-purple-600 rounded-2xl transform rotate-1"></div>
            <img
              src="https://images.unsplash.com/photo-1558618666-fcd25c85cd64?w=1200&h=600&fit=crop"
              alt="Smart home thermal management"
              className="relative rounded-2xl shadow-2xl w-full h-[400px] object-cover"
            />
          </div>
        </div>
      </section>

      {/* Features Section */}
      <section className="py-24 bg-white dark:bg-gray-900">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="text-center mb-16">
            <h2 className="text-4xl font-bold text-gray-900 dark:text-white mb-4">
              Intelligent Energy Optimization
            </h2>
            <p className="text-xl text-gray-600 dark:text-gray-400">
              Harness the power of AI to transform your energy management
            </p>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-blue-600 to-purple-600 rounded-full flex items-center justify-center mb-4">
                  <Brain className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  AI Analysis
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Advanced machine learning algorithms analyze your thermal data
                  and provide personalized optimization recommendations.
                </p>
              </div>
            </Card>

            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-green-600 to-teal-600 rounded-full flex items-center justify-center mb-4">
                  <TrendingDown className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  Cost Savings
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Reduce energy bills by up to 30% with intelligent thermal
                  management and real-time optimization suggestions.
                </p>
              </div>
            </Card>

            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-orange-600 to-red-600 rounded-full flex items-center justify-center mb-4">
                  <Thermometer className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  Smart Monitoring
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Track temperature, humidity, and power usage across all your
                  devices in real-time with intuitive dashboards.
                </p>
              </div>
            </Card>

            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-purple-600 to-pink-600 rounded-full flex items-center justify-center mb-4">
                  <BarChart3 className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  Analytics
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Comprehensive analytics and insights help you understand your
                  energy patterns and make informed decisions.
                </p>
              </div>
            </Card>

            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-indigo-600 to-blue-600 rounded-full flex items-center justify-center mb-4">
                  <Shield className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  Secure & Private
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Your data is encrypted and secure. We never share your
                  information with third parties.
                </p>
              </div>
            </Card>

            <Card>
              <div className="flex flex-col items-center text-center">
                <div className="w-16 h-16 bg-gradient-to-r from-yellow-600 to-orange-600 rounded-full flex items-center justify-center mb-4">
                  <Zap className="w-8 h-8 text-white" />
                </div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                  Instant Insights
                </h3>
                <p className="text-gray-600 dark:text-gray-400">
                  Get immediate feedback and recommendations. No waiting for
                  reports or manual analysis.
                </p>
              </div>
            </Card>
          </div>
        </div>
      </section>

      {/* CTA Section */}
      <section className="py-24 bg-gradient-to-r from-blue-600 to-purple-600">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <h2 className="text-4xl font-bold text-white mb-4">
            Ready to Optimize Your Energy?
          </h2>
          <p className="text-xl text-blue-100 mb-8">
            Join thousands of users saving money and reducing their carbon
            footprint with ThermalAI.
          </p>
          <Link to="/signup">
            <Button
              size="lg"
              className="bg-white text-blue-600 hover:bg-gray-100"
            >
              Start Free Trial
            </Button>
          </Link>
        </div>
      </section>

      {/* Footer */}
      <footer className="bg-gray-900 text-gray-400 py-12">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 text-center">
          <p>&copy; 2024 ThermalAI. All rights reserved.</p>
        </div>
      </footer>
    </div>
  );
}

```
