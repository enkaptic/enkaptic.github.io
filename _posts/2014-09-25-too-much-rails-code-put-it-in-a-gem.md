---
title: Too much Rails code? Put it in a gem!
tags:
- Computing
date: 2014-09-25 11:24:03 +0100
---
If you extract some of your Rails code into a gem, and test the gem separately, this will _guarantee_ that the gem cannot depend on the rest of your application. This is a Good Thing: if everything can depend on everything else, you end up with a hideous tangled mess.

And it's dead easy. First, assuming you are in your app's root directory, let's create the gem (using <code>-t</code> to include RSpec tests):

{% highlight bash %}
mkdir gems
cd gems
bundle gem my_gem -t
cd my_gem
{% endhighlight %}

First (of course), we're going to write a test and check that it fails. This goes in <code>spec/my_gem_spec.rb</code>. Here's the code:

{% highlight ruby %}
it 'should do something useful' do
  expect(Class.new.extend(MyGem).wibble).to eq("Wobble")
end
{% endhighlight %}

Check it fails:

{% highlight bash %}
rspec
# Failed examples:
#
# rspec ./spec/my_gem_spec.rb:8 # MyGem should do something useful
{% endhighlight %}

Then the code itself goes in <code>lib/my_gem.rb</code>. (You don't need to put all the code in there, of course; just put some <code>require</code> statements at the top, linking to the other files.) I've added a method to mine:

{% highlight ruby %}
module MyGem
  # Your code goes here...
  def wibble
    "Wobble"
  end
end
{% endhighlight %}

Check <code>rspec</code> works ... and that's it: you now have a gem!

The next task is to include <code>my_gem</code> in the Rails app. First add to the app's <code>Gemfile</code>:

{% highlight ruby %}
gem 'my_gem', path: 'gems'
{% endhighlight %}

Then run <code>bundle</code> in the app's root directory. (It will give some warning messages about the <code>gemspec</code>, but don't worry.) Then you can use the gem in your app. For example:

{% highlight ruby %}
class MyModel < ActiveRecord::Base
  include MyGem
end
{% endhighlight %}

Job done. Now <code>MyModel</code> has a <code>#wibble</code> method. Just what you needed:

{% highlight ruby %}
MyModel.new.wibble #=> "Wobble"
{% endhighlight %}
