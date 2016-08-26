---
title: WordPress to Rails (4) uploads
tags:
- Computing
redirect_to: http://www.anthonysmith.me.uk/2013/10/21/wordpress-to-rails-4-uploads/
---
WordPress allows you to upload files (such as images) and include them in posts. These end up in a location such as <code>.../wp-content/uploads/2012/06/githuboctocat.jpeg</code>.

I'm trying to keep things simple (for now), so I want all the uploaded files in one folder, with no fancy way of managing the files, and with only one version of each image (no thumbnails). So the above file would move to  <code>.../uploads/githuboctocat.jpeg</code>.

I didn't want this folder to be part of the application code, so (following <a href="http://stackoverflow.com/questions/1287342/is-it-safe-to-use-capistrano">this</a>) I added the following to <code>deploy.rb</code>:

{% highlight ruby %}
namespace :deploy do
  # ...
  after "deploy:setup", "deploy:uploads:setup"
  after "deploy:symlink", "deploy:uploads:symlink"

  namespace :uploads do
    desc "Create the uploads dir in shared path."
    task :setup do
      run "cd #{shared_path}; mkdir uploads"
    end

    desc "Link uploads from shared to current."
    task :symlink do
      run "cd #{current_path}/public; rm -rf uploads; ln -s #{shared_path}/uploads ."
    end
  end
end
{% endhighlight %}

Next, with the uploaded files in that new folder, the database needed to be migrated to change the URLs for the uploaded files. The following migration seems to work (and can be rolled back, as long as the uploaded files shares the same year and month as the post that refers to it):

{% highlight ruby %}
class ChangeWordpressUploadsUrls < ActiveRecord::Migration
  def up
    Post.all.each do |p|
      p.post_content = p.post_content.gsub(/http:\/\/www\.anthonysmith\.me\.uk\/research\/wp-content\/uploads\/\d{4}\/\d{2}\/([0-9a-zA-Z\-_.]+)/, '/uploads/')
      p.save
    end
  end

  def down
    Post.all.each do |p|
      # NB: assumes upload has same year and month as post
      year, month, day, name = p.permalink.split('/')
      p.post_content = p.post_content.gsub(/\/uploads\/([0-9a-zA-Z\-_.]+)/, "http://www.anthonysmith.me.uk/research/wp-content/uploads/#{year}/#{month}/" + '')
      p.save
    end
  end
end
{% endhighlight %}
